#!/usr/bin/env python
version = 1//1 ; ''' These two lines are
/*                   really important!  '''

# The first half of this program is written in python. It will count reads in BAM/SAM files using different metrics (per chromosome, per read length, etc).
# To learn more visit http://ac.gt/seqc or run "python ./SeQC.js.py --help"
# The second half of this program is written in JavaScript for using in Node.js. It will create a webserver that can interact with the database created by 
# the Python code. You can get more information via http://ac.gt/seqc, or run "node ./SeQC.js.py --help"


## One day i'll get around to putting this up on pypi and/or npm, but for now here's instructions on how to do things on OSX:

# Pysam Installation
########################
# For Mac OSX:
# > sudo easy_install pip
# > pip install --user pysam

# Postgres Installation:
########################
# For Mac OSX:
# > ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
# > brew install postgres
# > sudo easy_install pip
# > pip install --user psycopg2

## Import required librarys:
import os
import sys
import csv
import json
import time
import urllib
import getpass
import sqlite3
import hashlib
import StringIO
import argparse
import datetime
import fileinput
import subprocess
import collections

# Import optional librarys:
try: import pysam ; pysamInstalled = True
except: pysamInstalled = False
try: import psycopg2 ; postgresInstalled = True
except: postgresInstalled = False

## Import stat modules from SeQC's directory:
availableStats = {}
SeQC = os.path.abspath(__file__)                          ## Path to SeQC.js.py itself. Used for spawning new processes.
for potentialStat in os.listdir(os.path.dirname(SeQC)):
        thisFile = os.path.join(os.path.dirname(SeQC), potentialStat)
        if thisFile.endswith('.stat'):
            execfile(thisFile)                            ## I should upgrade this to a proper plugin/module manager, but this works fine for now. 

## Parse user-supplied command line options:
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
    description="Put in multiple BAM/SAM files, get out an SQL database of read statistics.")
parser.add_argument("-a", "--analysis", nargs='+', metavar='', action='append',
    help='''Optional. One or more statistic to gather.
Eg. "SeQC -a gc type -a tlen" would generate three stats: gc and type (linked) and tlen (unlinked).
Currently avalible stat names:
''' + '\n'.join(availableStats.keys()) + '''
Default is to do them all (linked)''')
parser.add_argument("--samtools", default='samtools', metavar='/path/to/samtools', 
    help="Optional until it isn't. Path to samtools if it's needed and cannot be found.")
parser.add_argument("-i", "--input", nargs='+', metavar='file',
    help='Required. One or more SAM (and/or BAM) files to analyse.')
parser.add_argument("-o", "--output", default='myProject.SeQC', metavar='',
    help='Optional. Name of output database (or path for SQLite). Default is "myProject".')
parser.add_argument("-q", '--quiet', action='store_true',
    help='Optional. No status bars in output. Good for logs, bad for humans.')
parser.add_argument("-pu", "--pguser", metavar='',
    help="Required if using Postgres. User account name.",)
parser.add_argument("-pp", "--pgpass", metavar='',
    help="Optional. Password of --pguser.")
parser.add_argument("-ph", "--pghost", default="localhost", metavar='',
    help="Optional. Hostname of postgres database. Default is localhost.")
parser.add_argument("--cpu", default=2, metavar='n', type=int,
    help="Optional. Number of processes/cores you want to use. Default is 2.")
parser.add_argument('--writeover', action='store_true',
    help="Optional. Will write over old data even if inputs/analyses are identical")
parser.add_argument("--debug", action='store_true',
    help="To err is human; to debug, divine.")
parser.add_argument("--SAM", help=argparse.SUPPRESS)                        # Used internally to tell subprocesses we are reading SAM. Set to either 'stdin' or 'file'.
parser.add_argument("--BAM", action="store_true", help=argparse.SUPPRESS)   # Used internally to tell subprocesses we reading directly via pysam or htspython.
args = parser.parse_args()

## Normalize user analysis arguments:
if args.analysis is not None:
    # We make each linked group of stats a set, before converting to a tuple, in the event that the user adds the same stat twice, eg. -a gc gc becomes just (gc).
    # We also sort by stat name in each linked group of stats so "-a tlen gc" becomes (gc,tlen).
    # We then put these tuples into a set, so that "-a gc tlen -a tlen gc" would become just ((gc,tlen)).
    # Finally we then put all these sorted tuples into a sorted list. ~ phew ~
    allAnalyses = set()
    for linkedGroup in args.analysis:
        linkedGroup = set(linkedGroup)
        for stat in linkedGroup:
            if stat not in availableStats.keys(): print '''
            ERROR: I do not know how to calculate the statistic "''' + stat + '''". Are you sure you typed it right? (case-sensitive)
            If the .stat file for this statistic is not in the same directory as SeQC, you may need to download it from http://ac.gt/seqc
            '''; exit()
        allAnalyses.add(tuple(sorted(linkedGroup)))
    args.analysis = sorted(allAnalyses)                         # A sorted list of sorted tuples
else:
    args.analysis = [tuple(sorted(availableStats.keys()))]      # Also a list, but of just 1 tuple - containing all the possible stats (sorted)

## This function checks to see if a fine is binary (BAM) or ASCII (SAM) by looking for a "null byte" (0x00) in the first 20Mb of the file.
def bamCheck(fileName):
    with open(fileName, 'rb') as xamfile:
        for i in range(10): # 10 tries to find a null byte in 2Mb chunks.
            section = xamfile.read(2048)
            if '\0' in section:            return True  # BAM file.
            if len(section) < 2048:
                if i == 0 and not section: return None  # Empty file.
                else:                      return False # SAM file.

## args.SAM and args.BAM are only ever present in subprocesses. The main parent python code (executed by the user) starts below:
if not args.SAM and not args.BAM:
    subprocesses = {}

    ## First we check for any updates to SeQC.js.py -- nothing is ever auto-installed, this is just a warning.
    try:
        latestVersion = urllib.urlopen("http://ac.gt/seqc/version").read()
        if latestVersion != str(version): print '''
---------------------------------------------------
A new version of SeQC is available!
Visit http://ac.gt/seqc for more details :)
---------------------------------------------------'''
    except: pass

    ## The message below is fired if the user does not provide any command line parameters at all when running SeQC:
    if not args.input: print '''
Hello :)
You must supply at least one SAM or BAM file, in the format "SeQC.js.py --input ./somefile1 ./somefile2 ... etc"
If this is your first time using SeQC, try "SeQC.js.py --help" or visit http://ac.gt/seqc/ for usage infomation.'''; exit()

    ## Determine if the user is writing to a SQLite or Postgres database, and determine that the provided details work:
    if args.pguser == None:
        print '   [ Using SQLite for output ]'
        if os.path.isdir(args.output): print '''
ERROR: You have provided an existing DIRECTORY for your output, but I need a file name!
Please use the exact name of the output file you want (you can always rename it later)'''; exit()

        if os.path.isfile(args.output): print '   [ Output file already exists ]'; newDB = False
        else:                           print '   [ Output file does not exist ]'; newDB = True
    else:
        print '   [ Using Postgres for output ]'
        if postgresInstalled:
            if args.pgpass == None:
                if args.debug: 
                    password = 'SOME_SECRET_PASSWORD' # we dont want --debug to print passwords! Since --debug on the parent never executes jobs anyway, this is fine.
                else:
                    if sys.stdin.isatty():
                        password = getpass.getpass('Enter password for account ' + args.pguser + ' : ')
                    else:
                        ## If you do not want to specify the password in the terminal with --pgpass for security reasons, but you want
                        ## to run SeQC as part of a pipeline with Postgres, you can cat/pipe a password from a chmod 700 file to SeQC,
                        ## and just " alias seqc='cat ./my/password.file | /path/to/SeQC.js.py --pguser yourAccount' "
                        password = sys.stdin.readline().rstrip()
            else: password = args.pgpass
        else: print '''
You have specified a postgres username (and therefore want to use postgres) but you have not installed psycopg2!
Please install it via "pip install --user psycopg2"'''; exit()

    ## Check that the files the user wants to analyse can be accessed:
    ## If a file becomes inaccesible later it will be skipped, but a check early on dosen't hurt...
    usedInputs = []
    for inputFile in args.input:
        if not os.path.isfile(inputFile) or not os.access(inputFile, os.R_OK): 
            print 'ERROR: The input file ' + str(inputFile) + ' could not be accessed. Are you sure it exists and we have permissions to read it? (continuing without it)'
        else: usedInputs.append(inputFile)

    ## Most parameters to all subprocesses are static and never change, such as the database connection details or --writeover. 
    ## Here we collect them all into 1 string called explicitStats:
    explicitStats = ''
    for group in args.analysis: explicitStats += ' --analysis ' + ' '.join(group)
    if args.writeover: explicitStats += ' --writeover'
    if args.pguser != None:
        explicitStats += ' --pguser "' + args.pguser + '" --pgpass "' + password + '" --pghost "' + args.pghost + '"'

    ## We only need samtools if we have one or more BAM files. Here we check all files with the bamCheck function, 
    ## and if a binary file is found, we check if we can find/execute samtools.
    DEVNULL = open(os.devnull, 'wb')
    if args.debug: STDERR = subprocess.STDOUT
    else: STDERR = DEVNULL
    for inFile in usedInputs:
        ##### THIS IS WHERE LOG.BIO SUPPORT FOR PRE-HASHING WOULD GO (GOLD SUPPORT - PLATINUM REQUIRES PROGRAM TO MD5 as it reads file).
        if bamCheck(inFile):
            if not args.samtools:
                exitcode = subprocess.call('samtools', stdout=DEVNULL, stderr=DEVNULL, shell=True)
                if exitcode != 1: # Samtools is weird, in that calling it with no parameters returns exitcode 1 rather than 0.
                    print '''
ERROR: You have tried to calculate statistics on a BAM file, but you have not provided 
a path to samtools (and "samtools" from the command line does not seem to work...)
Please specify a path to samtools with the --samtools parameter :)'''; exit()
                else: args.samtools = 'samtools'
            else:
                code = subprocess.call(args.samtools, stdout=DEVNULL, stderr=DEVNULL, shell=True)
                if code != 1: print '''
ERROR: You have tried to calculate statistics on a BAM file, but the path to samtools you have provided:
("' + args.samtools + '") does not work :(
Please check it and try again :)\n'''; exit()
            if not pysamInstalled: print 'INFO: You should install pysam if you can. It makes everything a lot faster!'
            break # As soon as we have 1 BAM file, we do the above checks - but only once.

    ## This function fires off the subprocesses which actually analyse the input BAM/SAM data.
    ## It will be called in a loop later as we process the input files.
    def doSeQC(inputFile):
        if bamCheck(inputFile):
            if pysamInstalled:
                subprocessCommand = 'python "' + SeQC + '" --input "' + inputFile + '" --BAM --output "' + args.output + '" ' + explicitStats
                if args.debug: print subprocessCommand # Dont run it, just print what would have been run.
                else: return subprocess.Popen(subprocessCommand, stdout=subprocess.PIPE, stderr=STDERR, shell=True, executable='/bin/bash')
            else:
                subprocessCommand = '"' + args.samtools + '" view "' + inputFile + '" | python "' +  SeQC + '" --input "' + inputFile + '" --output "' + args.output + '" --SAM stdin ' + explicitStats
                if args.debug: print subprocessCommand
                else: return subprocess.Popen(subprocessCommand, stdout=subprocess.PIPE, stderr=STDERR, shell=True, executable='/bin/bash')
        else:
            subprocessCommand = 'python "' + SeQC + '" --input "' + inputFile + '" --SAM file --output "' + args.output + '" ' + explicitStats
            if args.debug: print subprocessCommand
            else: return subprocess.Popen(subprocessCommand, stdout=subprocess.PIPE, stderr=STDERR, shell=True, executable='/bin/bash')

    ## The table schema for the INFO table. 
    ## The INFO table stores metadata on all the analysed files, as well as 
    ## which stats (or linked groups of stats) have been collected for each file.
    ## I am always open to more suggestion on metadata we should be collecting!
    infoTableCreation = '''CREATE TABLE "INFO" (
                            "sampleHash" TEXT PRIMARY KEY,
                            "analyses" TEXT,
                            "sampleFileName" TEXT, 
                            "samplePath" TEXT, 
                            "creationTime" TEXT, 
                            "sampleSize" TEXT, 
                            "analysisTime" TEXT, 
                            "header" TEXT 
                        )'''

    ## Create the INFO table (if not already present) for SQLite:
    if args.pguser == None:
        try:
            con = sqlite3.connect(args.output, timeout=120)
            cur = con.cursor()
            if newDB: print '   [ Database created ]'
            cur.execute("SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='INFO';")
            if cur.fetchone()[0] == 1:
                print '   [ Output file is a database and can be accessed! ]'
            else:
                print '   [ Creating INFO table.. ]'
                cur.execute(infoTableCreation)
                con.commit()
                print '   [ INFO table created! ]'
        except sqlite3.Error, e:
            print '''
ERROR: Something went wrong creating the output database. The exact error was:
''' + str(e) + '''
If you do not know what caused this error, please e-mail longinotto@immunbio.mpg.de with the error message and I will help you :)'''; exit()
        finally:
            if cur: cur.close()
            if con: con.close()
    ## Create the INFO table (if not already present) for Postgres:
    else:
        try:
            con = psycopg2.connect(dbname='postgres', user=args.pguser, host=args.hostname, password=password)   # We connect to the main (always there) postgres database, just to establish a connection. From there we can create our new database table.
            con.set_isolation_level(psycopg2.extensions.ISOLATION_LEVEL_AUTOCOMMIT)                              # Required when doing some actions like database creation.
            cur = con.cursor()
            cur.execute("SELECT 1 from pg_database WHERE datname='" + args.output + "';")                        # See if database already exists.
            if cur.fetchone() == None:                                                                           # If not, create it...
                print '   [ Creating database.. ]'
                cur.execute('CREATE DATABASE "' + args.output + '"')
                con.commit()
                print '\n   [ Database "' + args.output + '" created - adding INFO table.. ]'
                cur.close();con.close()
                con = psycopg2.connect(dbname=args.output,user=args.pguser,host=args.hostname,password=password) # ...and connect to it.
                cur = con.cursor()
                cur.execute(infoTableCreation)
                con.commit()
                print '   [ INFO table added! You\'re good to go! ]'
            else:
                cur.close();con.close()
                con = psycopg2.connect(dbname=args.output,user=args.pguser,host=args.hostname,password=password) # Connect to existing database.
                cur = con.cursor()
                cur.execute("SELECT * from information_schema.tables where table_name='INFO'")
                if cur.rowcount == 1:
                    print '   [ Database exists as does INFO table! You\'re good to go! ]'
                else:
                    print '   [ Database exists but INFO table does not! Creating INFO table.. ]'
                    cur.execute(infoTableCreation)
                    con.commit()
                    print '   [ INFO table added! You\'re good to go! ]'
            cur.close();con.close()
        except psycopg2.ProgrammingError, e:
            print '''
ERROR: Something went wrong creating the output database. The exact error was:
''' + str(e) + '''
If you do not know what caused this error, please e-mail longinotto@immunbio.mpg.de with the error message and I will help you :)'''; exit()
        finally:
            if cur: cur.close()
            if con: con.close()


    ## Fire off the initial subprocesses to start analysing the data!
    if not args.quiet:
        print '\nStarting initial processes to analyse input files (max ' + str(args.cpu)  + ' at a time):\n'
    if args.debug: print '''
-------------------------------------------------------------------------------------------------------------------------------
Because you are in debug mode, data wont actually be analysed. Below are just the subprocess commands that WOULD have been run.
You can run them manually, also with or without a "--debug", to see more information on the error you are experiencing.
-------------------------------------------------------------------------------------------------------------------------------
'''
    for index, inputFile in enumerate(usedInputs):
        if index == args.cpu: break
        subprocesses[inputFile] = doSeQC(usedInputs.pop(0))
        if not args.quiet and not args.debug:
            print '\033[A   [ ' + str(index+1) + '/' + str(args.cpu) + ' ]'
    if args.debug: exit()

    ## Function to tidy up finished/aborted subprocesses:
    def tidy(msg):
        if args.quiet: print msg
        else: allDone.append(msg)
        del subprocesses[theFile]
        try: del outputs[theFile]
        except KeyError: pass 

    ## Subprocesses always respond with 1 character of status data every 1 second (unless they have become unresponsive)
    ## This is because when I first wrote SeQC 2 years ago, I didnt know about select.select() and non-blocking stdout reading.
    ## In the future I will use select, but even though the below is ugly as hell, it's never failed me yet!
    ## Basically we draw the status bars, and submit new subprocesses when old ones complete.
    print '\n\nCalculating Statistics..'
    outputs = {}
    timeToWrite = []
    arrow = '>'
    allDone = []
    while len(subprocesses) > 0:
        ## Get the terminal width:
        terminalRows, terminalCols = map(int,os.popen('stty size', 'r').read().split())

        ## Before we print out the 'active' stuff, lets print out the completed/failed stuff so that it scrolls up and out of the way.
        for item in allDone: print item.ljust(terminalCols)  ## Using ljust here so we wipe over the stuff previously on this line.
        allDone = []

        ## Read the 1 byte of status output from all subprocesses:
        for theFile, handler in subprocesses.items():
            try: out = handler.stdout.read(1)
            except: tidy('ERROR: Process analysing ' + theFile + ' unexpectedly stopped?!')

            if out in ['0','1','2','3','4','5','6','7','8','9']: outputs[theFile][1] += int(out)
            elif out == '#': outputs[theFile] = ['Calculating MD5         |',0]                                     ## Statuses are always
            elif out == '|': outputs[theFile] = ['Calculating Statistics  |',0]                                     ## exactly 25 characters
            elif out == '@': outputs[theFile] = ['Writing To Database     |',0]                                     ## in width and include
            elif out == '$': outputs[theFile] = ['Postgres Table Indexing |',0]                                     ## the starting pipe "|"
            elif out == '%': tidy('File ' + theFile + ' skipped as it is already in the database.')
            elif out == '?': tidy('DATA ERROR: ' + theFile + ' aborted by SeQC. Rerun with --debug for more info.')
            elif out == '!': tidy('Completed: ' + theFile)
            else: tidy('SeQC ERROR: ' + theFile + ' failed, likely due to a bug in SeQC. Run with "--debug" for more info.')

        ## Print some beautiful status bars ;)
        progressBarSpace = int(terminalCols/2)                              ## Use half the screen for the progress bar and status.
        fileNameSpace = terminalCols - progressBarSpace                     ## Use the rest for the file name (full path).
        progressScaling = (progressBarSpace -27) / 100.                     ## -27 because we dont want the scaling factor to know about the status, arrow, or final pipe.
        arrow = '' if arrow else '>'                                        ## '' is falsey, while '>' is truthy, so this toggles it.
        for theFile, storedOut in outputs.items():
            if len(theFile) > fileNameSpace:
                left = '<-- ' + (theFile + ': ')[-fileNameSpace+4:]         ## <-- to indicate the filename was truncated
            else:
                left = theFile.ljust(fileNameSpace)
            right = (storedOut[0] + ('=' * int(storedOut[1] * progressScaling)) + arrow).ljust(progressBarSpace-1) + '|'
            if not args.quiet: print left + right

        if not args.quiet:
            for x in range(0,len(outputs)): sys.stdout.write('\033[A')  ## Put the cursor back up as many rows in the terminal as we have data files
        sys.stdout.flush()

        ## Run some more subprocesses if we need to:
        if len(subprocesses) < args.cpu:
            try:
                inputFile = usedInputs.pop(0);
                subprocesses[inputFile] = doSeQC(inputFile)
            except IndexError: pass ## No more files to process!

    ## Finished processing data!
    for item in allDone: print item.ljust(terminalCols) ## print out any outstanding messages.
    if not args.quiet:
        for x in range(0,len(outputs)+1): print ''  ## put the cursor at the bottom of the terminal.

    ## Due to the way SQLite is implimented, with only 1 process being able to write to the database at a time, indexing has to be done at the very end.
    ## Otherwize, the subprocesses will all have to wait while 1 process indexes, and this bottleneck slows things down considerably. It is not an issue
    ## for Postgres however, which uses multiple worker threads. Below we do the indexing for SQLite.
    if args.pguser == None:
        print 'All done calculating file statistics! Now creating indexes on your SQlite tables:'
        try:
            con = sqlite3.connect(args.output, timeout=120)
            cur = con.cursor()
            cur.execute("SELECT name FROM sqlite_master WHERE type='table';")   # get a list of all tables in our database
            rows = cur.fetchall()
            for row in rows:
                tableName = row[0]
                tableColumns = ','.join(tableName.split('_')[1:])
                if len(tableName) > 32:
                    print 'Checking index on ' + tableName + ' ...',
                    sys.stdout.flush()
                    cur.execute("SELECT COUNT(*) FROM sqlite_master WHERE type='index' AND name='" + tableName + "_INDEX';")
                    if cur.fetchone()[0] == 1:
                        print 'index already exists!'
                    else:
                        print ' creating index ... ',
                        cur.execute("CREATE INDEX IF NOT EXISTS \"" + tableName + "_INDEX\" on \"" + tableName + "\" (" + tableColumns + ",counts)")
                        print 'done!'
                else:
                    print 'Not indexing ' + tableName
        except sqlite3.Error, e:
            print '\nERROR: Something went wrong indexing the tables of the database. The exact error was:\n' + str(e) + '\nIf you do not know what caused this error, please e-mail longinotto@immunbio.mpg.de with the error message and I will help you :)\n'; exit()
        finally:
            if con: con.close()

    ## That's it! That's the parent process taken care of. The rest of the python code is for the subprocesses that actually collect
    ## the statistics from each file, and put them into the SQL database. Here we go...
    print '\nAll done! Have a lovely day :)\n'
    exit()

elif args.BAM or args.SAM:
    ## We are only ever here if we're running from a subprocess!
    inputFile = args.input[0]

    ## I wrote this little class to be a generic reporter class that you can
    ## just drop in to a loop and get back periodic % completion updates
    class pinger:
        def __init__(self,updateTime):
            self.total = None
            self.percent = None
            self.nextPing = None
            self.updateTime = updateTime
        def change(self,mode,total=1):
            self.total = float(total) # float so divisions return accurate floats
            self.percent = 0
            self.nextPing = self.updateTime + time.time()
            sys.stdout.write(mode); sys.stdout.flush()
        def pong(self,position=0):
            if time.time() > self.nextPing:
                percentDiff = int((position/self.total)*100) - self.percent
                if percentDiff > 0:
                    sys.stdout.write('9') if percentDiff > 9 else sys.stdout.write(str(percentDiff))
                    self.percent += percentDiff
                else: sys.stdout.write('0')
                sys.stdout.flush()
                self.nextPing = 1 + time.time()
    ping = pinger(1)

    ## Get both MD5 hash and read count at the same time. 
    def MD5andCount(inputFile,samtools):
        ## We want to know the readcount of the BAM/SAM file to make status bars, but running samtools view -c takes an unacceptibly long time. MD5 however is very fast.
        ## We can do both simultaniously, essentially requiring the same amount of time as running samtools -c on the file.
        ## But we can go one better, by only passing the first X byte to samtools -c, and then estimating the size of the full file from that sample.
        estimate = True
        sizeInBytes = os.path.getsize(inputFile)
        pipe = subprocess.PIPE
        ## samtools doesn't print the count if the file is truncated, so we're stuck with this garbage:
        if estimate: readCounting = subprocess.Popen('"'+ samtools +'" view -h - | "'+ samtools +'" view -c -',stderr=pipe,stdin=pipe,stdout=pipe,shell=True,executable='/bin/bash')
        else:        readCounting = subprocess.Popen('"'+ samtools +'" view -c -',                             stderr=pipe,stdin=pipe,stdout=pipe,shell=True,executable='/bin/bash')
        md5 = hashlib.md5()
        chunkSize = 128 * md5.block_size
        ping.change('#',sizeInBytes)
        with open(inputFile,'rb') as f: 
            for x,chunk in enumerate(iter(lambda: f.read(chunkSize), b'')): 
                md5.update(chunk)
                if not estimate: readCounting.stdin.write(chunk)
                elif x < 1000:   readCounting.stdin.write(chunk)
                ping.pong(x*chunkSize)
        readCounting.stdin.close()
        readCounting.wait()
        if estimate:
            sample = int(readCounting.stdout.read())
            bytesUsed = chunkSize * 1000.
            readCount = int( (sample/bytesUsed) * sizeInBytes * 1.07 )
        else:
            readCount = int(readCounting.stdout.read())
        return [ readCount , md5.hexdigest() ]
    totalReads, fileHash = MD5andCount(inputFile,args.samtools)

    ## Now we check if any analyses have been performed on this sample before.
    ## If they have and --writeover is not set, we dont do the analysis again.
    if args.pguser == None: con = sqlite3.connect(args.output, timeout=120)
    else:                   con = psycopg2.connect(dbname=args.output, user=args.pguser, host=args.pghost, password=args.pgpass)
    cur = con.cursor()
    ## We have two options to find out if the table exists - check the SQL database itself, or check the INFO table. Which is 'more correct'? I think checking the INFO table
    ## is, as if it's not there, it doesn't matter whats already in the database, we're expected to add it in over the top - but below is the code to check the database directly
    ## cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='" + fileHash + "';")                                                            ## SQLite
    ## cur.execute("SELECT 1 FROM pg_catalog.pg_class WHERE relkind = 'r' AND relname = '" + args.output + "' AND pg_catalog.pg_table_is_visible(oid) LIMIT 1") ## Postgres
    cur.execute("SELECT analyses, 'Robert\"); DROP TABLE students; --' FROM INFO WHERE \"sampleHash\"='" + fileHash + "';")
    result = cur.fetchone() # we return a "little bobby tables" because if you just ask for 1 column, SQLite returns a tuple but Postgres a string. With two fields they both return tuples.
    if result != None:
        existingAnalyses = json.loads(result[0])
        if not args.writeover:
            skip = set()
            for analysis in args.analysis:
                if analysis in [ tuple(existingAnalysis['stats']) for existingAnalysis in existingAnalyses ]: skip.add(analysis)
            for x in skip: args.analysis.remove(x)
            if len(args.analysis) == 0:
                sys.stdout.write('%'); sys.stdout.flush() # Skipping
                exit()
    else: existingAnalyses = []
    if cur: cur.close()
    if con: con.close()


    ## Finally, we now start analyzing the data...
    ping.change('|',totalReads)
    data = [collections.defaultdict(int) for x in range(0,len(args.analysis))] # a list of dictionaries, one for every analysis group.
    if args.SAM:
        for stat in availableStats:
            availableStats[stat].process = availableStats[stat].SAM
        if args.SAM == 'stdin':
            sys.argv = '' ## Done so fileinput takes stdin and not args.
            inputData = csv.reader(fileinput.input(), delimiter='\t')
        else: 
            inputData = csv.reader(fileinput.input(inputFile), delimiter='\t') # fileinput also can read a file on disk
        header = ''
    elif args.BAM:
        for stat in availableStats:
            availableStats[stat].process = availableStats[stat].BAM
        inputData = pysam.Samfile(inputFile, "rb")
        header = inputData.header

    for processedReads, line in enumerate(inputData):
        try:
            ping.pong(processedReads)
            for a, analysis in enumerate(args.analysis):
                thisLine = tuple([ availableStats[stat].process(line) for stat in analysis ])
                data[a][thisLine] += 1
        except IndexError:
            ## Extremely crude header-skipping for SAM files! You can't trust 'comment' symbols, because I've seen SAM files which don't use them.
            if len(line) < 11: header += line
            if args.debug: print 'Skipped line: ' + str(line)
            exit()

    ## All done reading file for stats.
    ping.change('@',len(args.analysis))

    if args.pguser == None:
        con = sqlite3.connect(args.output, timeout=120)
        con.isolation_level = 'EXCLUSIVE'
        con.execute('BEGIN EXCLUSIVE')
        cur = con.cursor()
    else:
        con = psycopg2.connect(dbname=args.output, user=args.pguser, host=args.pghost, password=args.pgpass)
        cur = con.cursor()

    analysisJSON = []
    for a, analysis in enumerate(args.analysis):
        # We start by building up the SQL "CREATE TABLE" and "INSERT" queries...
        listOfTuples = []
        pgMegaString = ''
        createColumn, insertColumn, columnQ = '','',''
        for stat in analysis:
            ping.pong(a)
            statType = availableStats[stat].SQL
            createColumn += (stat + ' ' + statType + ',')
            insertColumn += (stat + ',');
            if args.pguser == None: columnQ += '?,'
            else: columnQ += '%s,'

        ## Add the counts on the end of the above
        createColumn += 'counts INTEGER'
        insertColumn += 'counts'
        if args.pguser == None: columnQ += '?'   # use ? for SQLite
        else:                   columnQ += '%s'  # use %s for Postgres

        ## Check if we actually have any data (this is a real SAM/BAM file) before writing anything to the database:
        if len(data[a]) == 0:
            print 'ERROR: No data was collected from this file - if this is a binary file, it is not in BAM format. If it is a text file, it is not in SAM format!'
            exit()

        ## Chromsome is a pain, because it's not stored at the read-level per-se, rather it's converted from an int to the actual chromosome using a mapping in the header
        ## This means that if we were to bin reads based on their actual chromosome name, it would be really slow. Much faster to bin using the original ints, then convert
        ## to an actual chromosome name.... now.
        if 'chromosome' in analysis: chrIDX = analysis.index('chromosome')

        ## Iterate the data object:
        for thisTuple,value in data[a].items():
            ping.pong(a)
            ## poppin' tuples in the club
            tupleToList = list(thisTuple)
            tupleToList.append(value);

            ## Hastily fix broken things while no one is looking
            if args.BAM and 'chromosome' in analysis:
                if tupleToList[chrIDX] == -1: tupleToList[chrIDX] = 'Unmapped'
                else: tupleToList[chrIDX] = inputData.getrname(tupleToList[chrIDX])
            if args.SAM and 'chromosome' in analysis:
                if tupleToList[chrIDX] == '*': tupleToList[chrIDX] = 'Unmapped'

            ## If we are using SQLite:
            if args.pguser == None:
                returnOfTheTuple = tuple(tupleToList)
                listOfTuples.append(returnOfTheTuple)
            else:
                for index, thing in enumerate(tupleToList):
                    if thing == None: thing = '\N'      ## postgres identifies null as \N  -  note that its not a special symbol, it's literally just a backslash and a capital N. This can be changed with " WITH NULL AS 'yolo' " if you'd rather yolo was the signifier for null.
                    tupleToList[index] = str(thing)
                thisRow = '\t'.join(tupleToList)
                thisRow += '\n'
                pgMegaString += thisRow

        tableName = fileHash + '_' + '_'.join(analysis) 
        fileName = os.path.basename(args.input[0])
        filePath = os.path.abspath(args.input[0])[:-len(fileName)]
        creationTime = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(os.path.getmtime(args.input[0])))
        analysisTime = str(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')) ## Year-Month-Day Hour:Minute:Second
        sampleSize = os.path.getsize(args.input[0])
        tableSize = len(data[a])
        header = json.dumps(header)
        analysisJSON.append({
            'stats':analysis,
            'rows':tableSize,
            'types': [ availableStats[stat].SQL for stat in analysis ],
            'versions': [ availableStats[stat].VERSION for stat in analysis ],
            'parameters': [ availableStats[stat].PARAMETERS for stat in analysis ],
            'performed': analysisTime
        })
        del data[a]

        if args.pguser == None:
            cur.execute('DROP TABLE IF EXISTS "' + tableName + '"')
            cur.execute('CREATE TABLE "' + tableName + '"(' + createColumn + ')')
            cur.executemany('INSERT INTO "' + tableName + '"(' + insertColumn + ') values ( ' + columnQ + ' )', listOfTuples)
            con.commit()
        else:
            cur.execute('DROP TABLE IF EXISTS "' + tableName + '"')
            cur.execute("CREATE TABLE \"" + tableName + "\"(" + createColumn + ")")
            pgCopyTime = StringIO.StringIO(pgMegaString)
            cur.copy_from(pgCopyTime, '"' + tableName + '"')
            con.commit()


    '''
    OK John, just a tiny TINY bit left to go.
    Make sure pinger symbols and parent process match up.
    Can we know % of indexing done? (for SQLite and Postgres)
    A SERVER table in the database for telling node how to display the main data table page. Key/Value pairs.
    Needs to be flexible so people can have as many columns as they like and is involved in replies for info data table, but also somewhat standardized so that people can click on multiple
    rows of multiple front-ends to compare data.
    '''

    ## INFO TABLE ADDITION / MODIFICATION
    if len(existingAnalyses) == 0:
        analysisJSON = json.dumps(analysisJSON)
        if args.pguser == None: placeholder = '(?, ?, ?, ?, ?, ?, ?, ? )'
        else:                   placeholder = '(%s,%s,%s,%s,%s,%s,%s,%s)'
        cur.execute('INSERT INTO "INFO" ("sampleHash","analyses",    "sampleFileName","samplePath","creationTime","sampleSize","analysisTime", "header") VALUES ' + placeholder,
                                        ( fileHash,    analysisJSON,  fileName,        filePath,    creationTime,  sampleSize,  analysisTime,   header )                        )
    else:
        ## An INFO row already exists for this sample...
        for existingAnalysis in existingAnalyses:
            if tuple(existingAnalysis['stats']) not in args.analysis: analysisJSON.append(existingAnalysis)
        analysisJSON = json.dumps(analysisJSON)
        if args.pguser == None: sql = 'UPDATE "INFO" SET "analyses"=?, "sampleFileName"=?, "samplePath"=?, "analysisTime"=?, "header"=? WHERE "sampleHash"=?'
        else:                   sql = 'UPDATE "INFO" SET "analyses"=%s, "sampleFileName"=%s, "samplePath"=%s, "analysisTime"=%s, "header"=%s WHERE "sampleHash"=%s'
        cur.execute(sql,                                 (analysisJSON,  fileName,            filePath,        analysisTime,      header,           fileHash ))
    # cur.execute('DROP INDEX IF EXISTS "' + tableName + '_INDEX"')
    # cur.execute('CREATE INDEX "' + tableName + '_INDEX" on "' + tableName + '" (' + ','.join(linkedStats) + ',counts )')
    con.commit()
    cur.close()
    con.close()
    sys.stdout.write('!'); sys.stdout.flush()
 
        
'''*/

// I found some kind of weird bug in Node that if you type around 255 characters into the terminal, you fill up the STDIN pipe and stop receiving STDOUT. This is kind of dumb, so the two blocks below fix this:
process.stdin.on('readable', function() {
  var chunk = process.stdin.read();
  if (chunk !== null) {
    process.stdout.write('');
  }
});
process.stdin.on('end', function() {
    process.stdout.write('end');
});

// Set our global variables, required modules, etc :
try {
    fs = require('fs');
    async = require('async');
    http = require('http');
    argv = require('minimist')(process.argv.slice(2)); // slice 2 because the first 2 arguments are 'node' and 'SeQC.js.py'
    express = require('express');
    compress = require('compression');
    restapi = express();
    crypto = require('crypto');
    sqlite3 = require('sqlite3').verbose();
    options = { host: 'ac.gt', port: 80 }; // Main ac.gt website where you download the most up-to-date vizulization code.
    SeQCPath = '/seqc/' + version + '/'; // This is the path on ac.gt to grab your viz code from
    pg = require('pg');
}
catch (e) {
    console.log('ERROR: You are missing one or more of the required modules.Please run: \n"npm install async minimist express compression sqlite3 pg"\nto get the required modules. If that does not work, try:\n"npm --user install async minimist express compression sqlite3 pg" :) ')
    console.log(e)
    process.exit()
}
// It was originally decided that SeQC would maintain a persistent cache of query results (which is why we are so particular about normalizing the sql query) - but although it worked, complexity blew up out of control. To keep future bugs to a minumum, i removed it all (since its fast enough in my opinion anyway..) - but it might make a comeback as a non-persistent cache (just an object in Node which is dropped when the instance is terminated) - but even that isnt as trivial as you might think. We have to manage the size of this cache (its VERY hard to find out how much RAM a specific JS object uses...), and then we have a problem of different Node instances run in a cluster all having their own (not-shared) cache. But yeah, it might make a comeback in the next version of SeQC :)
// cacheSQL =           "CREATE TABLE IF NOT EXISTS \"CACHE\"     (querySampleHash TEXT, sampleHash TEXT, result TEXT, UNIQUE(querySampleHash) )";
// sessionSQL =         "CREATE TABLE IF NOT EXISTS \"SESSIONS\"      (id TEXT, samples TEXT, UNIQUE(id) )";
// sessionDataSQL =     "CREATE TABLE IF NOT EXISTS \"SESSIONDATA\" (id TEXT, analysis TEXT, lookingAt TEXT, sessionHashes TEXT, sessionPlots TEXT, UNIQUE(id) )";
cachedQueries = {};

// Get user command-line arguments:
if (argv['_'] == []) { console.log('\nHello - and well done for getting Node.js installed and working!\n\nTo run SeQC on your local computer, run:\n  "node SeQC.js.py /path/to/database.SeQC"\n\nIf you want others on your network to be able to access your SeQC database, you can change the hostname and port of the service by running:\n  "node SeQC.js.py /path/to/database.SeQC --port 8080 --hostname computer.name.or.ip"\n\n Best of Luck! :)\n'); process.exit(); }
else { databaseFile = argv['_'][0]; }               // The '_' array is where minimist stores arguments which didnt have a -- before them. For us, thats all the linked databases the user specified. Future versions of SeQC will allow multiple databases to be connected to your webservice - including remote ones! (like a publicly hosted ENCODE QC database hosted on ac.gt, etc)

if (argv.hostname) { hostname = argv.hostname; }    // --h or --hostname to specify hostname of your webservice
else if (argv.h) { hostname = argv.h; }
else { hostname = 'localhost'; }

if (argv.port) { port = argv.port; }                // --p or --port to specify the port of your webservice
else if (argv.p) { port = argv.p; }
else { port = '8080'; }

if (argv.pguser) { pguser = argv.pguser; }  // --pu or --pguser to specify the username to connect to your postgres database
else if (argv.pu) { pguser = argv.pu; }
else { pguser = 'postgres'; }

if (argv.pgpass) { pgpass = argv.pgpass; }  // --pp or --pgpass to specify the password to connect to your postgres database
else if (argv.pp) { pgpass = argv.pp; }
else { pgpass = ''; }

if (argv.pghost) { pghost = argv.pghost; }  // --ph or --pghost to specify the hostname of your postgres database
else if (argv.ph) { pghost = argv.ph; }
else { pghost = 'localhost'; }

if (argv.cpu) { var cpu = argv.cpu; }               // --cpu to specify how many postgres workers you want to use in your pool. This number should be equal to the maximum number of samples you wish to analyse at once. Here is some very good advice from Brian Carlson, author of node-postgres: 
else { var cpu = 10; }                              // " A connection to postgres takes a non-trivial amount of RAM on the server. I'm not sure the exact amount but somewhere around 1,000 connections on our prod instance we run out of memory completely and the box stops responding.  So...we try to always keep connections under 500 just to be safe and put an AWS alarm on the box if it goes above 600."
                                                    // " So you can tune your pool to be quite a bit larger. Another thing you can do if you have a 40 core machine is use a module like http://github.com/brianc/node-forky to fork your node process onto all 40 cores.  This will give you a almost fully linear 40x speedup BUT one thing to be careful about is you will have 40 processes. "
                                                    // " Each process will have its own connection pool. So if your poolSize is 20 you'll be using 800 connections which can be dangerous. You can tweak the pool size and concurrency level of forky to get some serious throughput though. "

exists = fs.existsSync(databaseFile);
if (exists) { 
    console.log('\nFound your database file! Attempting to start the webservice...');
    db = new sqlite3.cached.Database(databaseFile);
    db.serialize(function() {                                                   // Do i need to serialize this?
        db.run("PRAGMA case_sensitive_like = true",function(err) {
            if (err) { return initialSQLError(err); }
            else { console.log('... we do! Now just go to\n\n http://' + hostname + ':' + port + '/ \n\nin your browser to start analysing your data! :)\n'); postgresConnection = false; }
        });
    });
} else { 
    postgresConnection = 'postgres://' + pguser + ':' + pgpass + '@' + pghost + '/' + databaseFile
    pg.connect(postgresConnection, function(err, client, done) {
        if(err) {
            console.log(err)
            console.log('\nWHOOPS: The database you have specified does not exist! Please double-check the name or path and try again :)\n');
            process.exit();
        } else { 
            console.log('... we do! Now just go to\n\n http://' + hostname + ':' + port + '/ \n\nin your browser to start analysing your data! :)\n');
        }

        done();
    });
}

// Function to add objects to an existing object of objects, making branches in the tree if required.
objectify = function (o, l, v){
    var pl = l.pop(); 
    l.reduce(function(o,k){ 
        if(!(k in o)) o[k] = {};
        return o[k] 
    },o)[pl] = v;
    return o 
};

// Got SQL errors? Theres an app for that!
initialSQLError = function(err) {
    console.log('\nWHOOPS: While it appears you DID specify a file to use as your database, and this file does exist, either we cannot access it due to permissions, or it was never actually a database file in the first place!\nFull error: (' + err + ')\n')
    process.exit()
}
SQLError = function (error, rows) {
    if (error) {
        console.log("Error in callback");
        console.log(error);
    }
}

// Puts together JSON from the user's POST request as it comes in chunk by chunk
parseJSON = function(req, callback) {
    if (req.method == 'POST') {
        var jsonString = '';
        req.on('data', function (data) { jsonString += data; });
        req.on('end', function () {
            callback(JSON.parse(jsonString));
        });
    } else { res.json('fail.'); }
}

// Case insensitive sorter used to sort read flag filters
caseInsensitiveSort = function(x,y) {
  if (x.toLowerCase() > y.toLowerCase()){ return 1; }
  return 0;
};

// Takes random flag filters in any order and sorts/groups them to make the most efficient SQL query. So 'KcB' becomes 'Bc' and 'K' (as you'll see in the SQL commands used) :
optimizeFlagFilters = function (letters) {
  var output   = [];
  var lastMemo = letters.reduce(function(memo, letter) {
    var lastIndex      = memo.length - 1;
    var nextCharCode   = memo.toLowerCase().charCodeAt(lastIndex) + 1;
    var letterCharCode = letter.toLowerCase().charCodeAt(0);
    if (nextCharCode === letterCharCode) { return memo + letter; }
    output.push(memo);
    return letter;
  });
  output.push(lastMemo);
  return output;
}

// Takes query conditions via user's JSON POST and returns a list sanitized/normalized SQL queries :
sanitizeQuery = function (queryConditions) {
    var userData = []; 
    // Normalize/sanitize SELECT and GROUP:
    if ( ('subplotOn' in queryConditions) || ('lookingAt' in queryConditions) ) {
        if (queryConditions.subplotOn == undefined)         { return 'fail' }
        if (queryConditions.subplotOn == "chromosome")      { sqlSELECT = "SELECT sum(counts),chromosome"   ; sqlGROUP = " GROUP BY chromosome" }
        if (queryConditions.subplotOn == "tlen")            { sqlSELECT = "SELECT sum(counts),tlen"         ; sqlGROUP = " GROUP BY tlen"       }
        if (queryConditions.subplotOn == "type")            { sqlSELECT = "SELECT sum(counts),type"         ; sqlGROUP = " GROUP BY type"       }
        if (queryConditions.subplotOn == "gc")              { sqlSELECT = "SELECT sum(counts),gc"           ; sqlGROUP = " GROUP BY gc"         }
        if (queryConditions.subplotOn == "flag")            { sqlSELECT = "SELECT sum(counts),flag"         ; sqlGROUP = " GROUP BY flag"       }
        if (queryConditions.subplotOn == "None")            { sqlSELECT = "SELECT sum(counts)"              ; sqlGROUP = ""                     }

        if (queryConditions.lookingAt == undefined)     { return 'fail' }
        if (queryConditions.lookingAt == "chromosome")  { sqlSELECT += ",chromosome FROM\"" ; if (queryConditions.subplotOn == "None"){ sqlGROUP = " GROUP BY chromosome ;"     } else { sqlGROUP += ",chromosome ;"    };  };
        if (queryConditions.lookingAt == "tlen")        { sqlSELECT += ",tlen FROM \""      ; if (queryConditions.subplotOn == "None"){ sqlGROUP = " GROUP BY tlen ;"           } else { sqlGROUP += ",tlen ;"          };  };
        if (queryConditions.lookingAt == "type")        { sqlSELECT += ",type FROM \""      ; if (queryConditions.subplotOn == "None"){ sqlGROUP = " GROUP BY type ;"           } else { sqlGROUP += ",type ;"          };  };
        if (queryConditions.lookingAt == "gc")          { sqlSELECT += ",gc FROM \""        ; if (queryConditions.subplotOn == "None"){ sqlGROUP = " GROUP BY gc ;"             } else { sqlGROUP += ",gc ;"            };  };
        if (queryConditions.lookingAt == "flag")        { sqlSELECT += ",flag FROM \""      ; if (queryConditions.subplotOn == "None"){ sqlGROUP = " GROUP BY flag ;"           } else { sqlGROUP += ",flag ;"          };  };
        if (queryConditions.lookingAt == "counts")      { sqlSELECT += " FROM \""           ; sqlGROUP += " ;"                                                                                                              };
    } else { return 'fail' } ;


    // Normalize/sanitize WHERE:
    var allFilters = [];

    if ('chromosome' in queryConditions) {
        // check chromosome syntax
        var chromosomeFilters = "chromosome in ('";
        var chromosomes = queryConditions.chromosome;
        for (var i = 0; i < chromosomes.length; i++) { chromosomes[i] = String(chromosomes[i]); }
        chromosomes.sort();
        chromosomeFilters += chromosomes.join("','") + "')"
        allFilters.push(chromosomeFilters);
    }

    if ('type' in queryConditions) {
        // check type syntax
        var typeINFilters = [];
        var type = queryConditions.type;
        for (var i = 0; i < type.length; i++) {
            var filter = String(type[i]).split(':');
            if (filter.length == 1) {
                if (!isNaN(filter[0]) && (filter[0] >= 0) && (filter[0] <= 20) ) { 
                    typeINFilters.push(filter[0]); 
                } else { return 'fail' } // Not a number between 0 and 20
            } else { return 'fail' } // Cant/wont do ranges!
        };

        if (typeINFilters.length == 1) { allFilters.push('type=' + typeINFilters[0]); }
        else if (typeINFilters.length > 1) {
            typeINFilters.sort();
            allFilters.push('type in (' + typeINFilters.join(',') + ')' );
        }
    }

    if ('tlen' in queryConditions) {
        var tlenINFilters = [];
        var tlenBETWEENFilters = [];
        var tlenBOTHFilters = [];
        var tlen = queryConditions.tlen;
        for (var i = 0; i < tlen.length; i++) {
            var filter = String(tlen[i]).split(':'); // if you try and .split() a number, you are going to have a bad time..
            if (filter.length == 1) {
                if (!isNaN(filter[0])) { tlenINFilters.push(filter[0]); } 
                else { return 'fail' } // Not a number
            } else if (filter.length == 2) {
                if ( !isNaN(filter[0]) && !isNaN(filter[1]) ) {
                    if ( Number(filter[0]) < Number(filter[1]) ) { tlenBETWEENFilters.push('(tlen BETWEEN ' + filter[0] + ' AND ' + filter[1] + ')'); }
                    else if ( Number(filter[0]) == Number(filter[1]) ) { tlenINFilters.push(filter[0]); } // both numbers the same.
                    else { return 'fail' } // FROM bigger than TO
                } else if ( (filter[0] == 'min') && Number(filter[1]) ) {
                    tlenBETWEENFilters.push('tlen < ' + filter[1]);
                } else if ( Number(filter[0]) && (filter[1] == 'max') ) {
                    tlenBETWEENFilters.push('tlen > ' + filter[0])
                } else { return 'fail' } // Not two numbers
            } else { return 'fail' } // something like a:b:c
        };

        if (tlenINFilters.length == 1) { tlenBOTHFilters.push('tlen=' + tlenINFilters[0]); }
        else if (tlenINFilters.length > 1) {
            tlenINFilters.sort();
            tlenBOTHFilters.push('tlen in (' + tlenINFilters.join(',') + ')' );
        }

        if (tlenBETWEENFilters.length >= 1) { 
            tlenBOTHFilters.push(tlenBETWEENFilters.join(' OR '));
        }
        allFilters.push(tlenBOTHFilters.join(' OR '));
    }

    if ('gc' in queryConditions) {
        // check gc syntax
        var gcINFilters = [];
        var gcBETWEENFilters = [];
        var gcBOTHFilters = [];
        var gc = queryConditions.gc;
        for (var i = 0; i < gc.length; i++) {
            var filter = String(gc[i]).split(':');
            if (filter.length == 1) {
                if (!isNaN(filter[0]) && (filter[0] >= 0) && (filter[0] <= 100) ) { 
                    gcINFilters.push(filter[0]); 
                } else { return 'fail' } // Not a number between 0 and 100
            } else if (filter.length == 2) {
                if ( !isNaN(filter[0]) && !isNaN(filter[1]) && (filter[0] >= 0) && (filter[0] <= 100) && (filter[1] >= 0) && (filter[1] <= 100) ) {
                    if ( Number(filter[0]) < Number(filter[1]) ) { gcBETWEENFilters.push('(gc BETWEEN ' + filter[0] + ' AND ' + filter[1] + ')'); }
                    else if ( Number(filter[0]) == Number(filter[1]) ) { gcINFilters.push(filter[0]); } // both numbers the same...
                    else { return 'fail' } // FROM bigger than TO
                } else { return 'fail' } // Not two numbers between 0 and 100
            } else { return 'fail' } // something like a:b:c
        };

        if (gcINFilters.length == 1) { gcBOTHFilters.push('gc=' + gcINFilters[0]); }
        else if (gcINFilters.length > 1) {
            gcINFilters.sort();
            gcBOTHFilters.push('gc in (' + gcINFilters.join(',') + ')' );
        }
        if (gcBETWEENFilters.length >= 1) { 
            gcBOTHFilters.push(gcBETWEENFilters.join(' OR '));
        }
        allFilters.push(gcBOTHFilters.join(' OR '));
    }

    if ('flag' in queryConditions) {
        // check flag syntax
        var flagFilters = [];
        var flag = queryConditions.flag;
        for (var i = 0; i < flag.length; i++) {
            if ( /^[a-lA-L]{1,12}$/.test(flag[i]) ) {
                var optimizedFlags = flag[i].split('').sort(caseInsensitiveSort);
                optimizedFlags = optimizeFlagFilters(optimizedFlags);
                for (var f = 0; f < optimizedFlags.length; f++) {
                    optimizedFlags[f] = "flag LIKE '%" + optimizedFlags[f] + "%'";
                };
                flagFilters.push(optimizedFlags.join(' AND '))
            } else { return 'fail' } // Not flags
        }
        allFilters.push('(' + flagFilters.join(') OR (') + ')');
    }
    var sqlWHERE = '\" WHERE (' + allFilters.join(') AND (') + ')';
    if (sqlWHERE == '\" WHERE ()') { sqlWHERE = '\"'; }

    return {SELECT: sqlSELECT, WHERE: sqlWHERE, GROUP: sqlGROUP};
}

// This is where the magic happens
doQuery = function(queryConditions, sql, res, startTime) {


    if (sql.WHERE == '\"') { queryHash = 'Unfiltered' }                                                                             // queries with no filters are useful for normalizing things on the webpage :)
    else { queryHash = crypto.createHash('md5').update([sql.SELECT,sql.WHERE,sql.GROUP].join('')).digest("hex").substring(0,4); }   // queryHash is specific to the query, irrispective of the samples analysed. Also its only 4 characters long because, lets face it, no ones going to do a billion unique queries.

    if (queryConditions.subplotOn == 'tlen' || queryConditions.lookingAt == 'tlen') {
        var dumpTlenZero = true
    } else {
        var dumpTlenZero = false
    }

    var results = {};
    var startQueryTimes = {};
    var endQueryTimes = {};
    async.each(queryConditions.samples, function(sample,callback){
        var sampleQuery = [sql.SELECT,sample,sql.WHERE,sql.GROUP].join(''); // This is the actual query!
        var querySampleHash = crypto.createHash('md5').update(sampleQuery).digest("hex").substring(0,6); // OK, no one is going to do a billion unique sample-queries either, but having 4-character strings for queries and 6-character strings for samples-queries makes it easier to identify what a hash is as soon as you look at it :)
        startQueryTimes[sample] = [sampleQuery, new Date().getTime()];
        if (postgresConnection) {
            pg.connect(postgresConnection, function(err, client, done) {
                if(err) {
                    return console.error('Error fetching client from pool', err);
                }
                client.query(sampleQuery, function(err, result) {
                    done();
                    if(err) {
                        return console.error('error running query', err);
                    }
                    result.rows.forEach(function (row,key,y) {
                        if (row[queryConditions.subplotOn] == undefined) { var analysisCatagory = 'None' } else { var analysisCatagory = row[queryConditions.subplotOn] }
                        if (row[queryConditions.lookingAt] == undefined) { var lookingAtCatagory = 'counts' } else { var lookingAtCatagory = row[queryConditions.lookingAt] }
                        if (dumpTlenZero && (row.tlen === 0 || row.sum < 10) ) { return }   // The presence of a tlen value of 0 when there are no actual fragments of length 0 (as is often the case as aligners use 0 as a way to say 'no fragment length') will screw up all the % transormations used in the subsequent plots. As much as it pains me to 'hide' data from scientists, more harm will come from keeping the 0 tlens in that removing them. How could it possibly be that in a format for storing data the number '0' stands for something other than the number zero? Well, its because computer scientists write aligner software, not biologists. 99% of the time this is a very good thing - but in some cases they leave their footprints (like 0 being an O.K. way to represent false) all over the data. In future versions of SeQC this will be fixed by checking, as the SAM format perscribes, that the read and mate flags are mapped, etc - but then we move further and further away from raw data analysis and into the realms of processed data analysis - along with all of my personal "unknown assumptions" about what that processing should involve - making SeQC less of a QC tool for finding 'bugs' in your sample prep and/or alignment, and more of a general Analysis tool - anyway, I think thats quite enough philosophy for one comment...
                        objectify(results, [analysisCatagory,queryHash,sample,lookingAtCatagory] , Number(row.sum) );
                    });
                    endQueryTimes[sample] = new Date().getTime();
                    callback();
                });
            });
        } else {
            // HERE IS WHERE YOU TEST SQLITE PARALLEL CONNECT - init db.
            db.all(sampleQuery, function(err,rows){
                rows.forEach(function (row,key,y) {
                    if (row[queryConditions.subplotOn] == undefined) { var analysisCatagory = 'None' } else { var analysisCatagory = row[queryConditions.subplotOn] }
                    if (row[queryConditions.lookingAt] == undefined) { var lookingAtCatagory = 'counts' } else { var lookingAtCatagory = row[queryConditions.lookingAt] }
                    if (dumpTlenZero && (row['tlen'] === 0 || row['sum(counts)'] < 10) ) { return } // The presence of a tlen value of 0 when there are no actual fragments of length 0 (as is often the case as aligners use 0 as a way to say 'no fragment length') will screw up all the % transormations used in the subsequent plots. As much as it pains me to 'hide' data from scientists, more harm will come from keeping the 0 tlens in that removing them. How could it possibly be that in a format for storing data the number '0' stands for something other than the number zero? Well, its because computer scientists write aligner software, not biologists. 99% of the time this is a very good thing - but in some cases they leave their footprints (like 0 being an O.K. way to represent false) all over the data. In future versions of SeQC this will be fixed by checking, as the SAM format perscribes, that the read and mate flags are mapped, etc - but then we move further and further away from raw data analysis and into the realms of processed data analysis - along with all of my personal "unknown assumptions" about what that processing should involve - making SeQC less of a QC tool for finding 'bugs' in your sample prep and/or alignment, and more of a general Analysis tool - anyway, I think thats quite enough philosophy for one comment...
                    objectify(results, [analysisCatagory,queryHash,sample,lookingAtCatagory] , row['sum(counts)'] );
                });
                endQueryTimes[sample] = new Date().getTime();
                callback();
            });
        }

    }, function(err){
        var endDoAllQueries = new Date().getTime();
        var timeDoAllQueries = endDoAllQueries - startTime;
        res.json(results);
        var endRequest = new Date().getTime();
        var timeDoRequest = endRequest - startTime;
        console.log('Request for SQL data:')
        Object.keys(startQueryTimes).forEach(function(sample) { 
            var thisTime = endQueryTimes[sample] - startQueryTimes[sample][1]
            console.log('"' + startQueryTimes[sample][0] + '"  -  ' + thisTime + 'ms')
        });
        console.log('Done talking to database in ' + timeDoAllQueries + 'ms')
        console.log('Whole request took ' + timeDoRequest + 'ms');
    });
};

// Actual Webserver Settings:

/*
Basically, the webservice is a very simple, lightweight webserver that understands only 3 queries (yourhostname:8080/somepage.html is a query in the world of webservers)
These queries are:

/samples        - Returns a list of all samples in your database (from INFO table) to populate the samples table on the homepage.
/updateSamples  - Updates a sample name or project name by writing to the INFO table.
/getData        - Returns the result of a specific SQL query (in JSON format) for use in vizulizations.

If the webservice gets a query other than one of the above, it will try two things:

First - It will look in the folder that you ran the webservice from (the folder you were in when you typed "node SeQC.js.py ..." ) and return a file if it sees one.
      - It will not process the file like a regular webserver (for example, it will not run PHP, etc), so you do not need to worry about code execution beyond what is in this file.
      - It will not accept queries which are above the directory it was run from, for example "webservice:8080/../../some_secret_file.txt" would NOT be allowed!
      - It will not return files which the user who ran node did not have access to - so you can run node from an underprivilaged account if you prefer.

Second - If no file was found on the local file system, it will query http://ac.gt/[version]/[query] for the file. If it gets a reply, it will pass it on to the user.
       - Absolutely no data is sent to ac.gt at any time. This is absolutely critical because many users of SeQC use it to analyse human genomic data, for which strict data protection laws exist across the world.

Why all this complexity? Well actually this makes updating/improving the SeQC project remarkably simple:
- We do not need to package a whole website into SeQC.js.py, which would otherwise need to be unpacked to a specific directory, permissions corrected, etc etc.
- The vizulization code gets updated VERY frequently as it is by far the most complicated code in the project. As a user of SeQC, you do not need to do anything to get the latest code! Just refresh the webpage :)
- SeQC takes your data privacy incredibly seriously - because the above three queries are the only queries an attacker could possibly use to talk to your server, the codebase is small and the risk of code execution is zilch. There are no mammoth PHP/CGI scripts which would take a week to deduce what they are actually doing. All the code executed on the server is in this file and its all very very simple!
- To run an 'offline' instance of SeQC, you just have to copy the code from http://ac.gt/seqc/1/ to your local directory, and your webservice will use that without having to visit ac.gt first
- Writing your own custom vizulization code is easy! Just copy the template you want to edit from http://ac.gt/seqc/1/ (say, SeQC.css), put it in the local node folder, and edit it! Dont like the colours? Choose your own! 
- If you want to share your custom vizulizations with the world, e-mail john (longinotto@immunbio.mpg.de) with your javascript/css code, and I will host it on http://ac.gt/seqc/ on a different visualization version number. To access it, just change your visualization version number at the very top of this file!

*/

restapi.use(express.static(__dirname, { maxAge:1, expires:1 })); // If you are frequently updating local code and do not want to cache anything, change maxAge and expires to 1ms :)
restapi.use(compress());  

restapi.get('/samples', function(req, res){
    var sampleRequestStart = new Date().getTime();
    if (postgresConnection) {
        pg.connect(postgresConnection, function(err, client, done) {
            if(err) {
                return console.error('Error fetching client from pool', err);
            }
            client.query('SELECT * FROM "INFO"', function(err, result) {
                done();
                if(err) {
                    return console.error('error running query', err);
                }

                allSamples = result.rows.map(function (row,key,y) { return { 'sampleHash':row.sampleHash, 'sampleID':row.sampleID, 'sampleFileName':row.sampleFileName, 'projectName':row.projectName, 'samplePath':row.samplePath, 'creationTime':row.creationTime, 'sampleSize':row.sampleSize, 'analysisTime':row.analysisTime, 'tableSize':row.tableSize, 'colour':row.colour }; });
                res.json(allSamples);
                var sampleRequestDuration = new Date().getTime() - sampleRequestStart;
                console.log('Sample List Requested - Time taken to reply: ' + sampleRequestDuration + 'ms');
            });
        });
    } else {
        db.all('SELECT * FROM "INFO"', function(err, rows){
            allSamples = rows.map(function (row,key,y) { return { 'sampleHash':row.sampleHash, 'sampleID':row.sampleID, 'sampleFileName':row.sampleFileName, 'projectName':row.projectName, 'samplePath':row.samplePath, 'creationTime':row.creationTime, 'sampleSize':row.sampleSize, 'analysisTime':row.analysisTime, 'tableSize':row.tableSize,  'colour':row.colour }; });
            res.json(allSamples);
            var sampleRequestDuration = new Date().getTime() - sampleRequestStart;
            console.log('Sample List Requested - Time taken to reply: ' + sampleRequestDuration + 'ms');
        });
    }
});

restapi.get('/updateSamples', function(req, res){
    if ( /^[a-f0-9]{32}$/.test(req.query['hash'])  && ( (req.query['column'] == 'sampleID') || (req.query['column'] == 'projectName') ) && /^[a-zA-Z0-9 _.():-]{1,50}$/.test(req.query['newVal']) ) {
        values = [ req.query['newVal'], req.query['hash'] ];
        console.log( 'UPDATE "INFO" SET "' + req.query['column'] + '"=\'' + String(values[0]) + '\' WHERE "sampleHash"=\'' + String(values[1]) + "'" )
        if (postgresConnection) {
            pg.connect(postgresConnection, function(err, client, done) {
                if(err) {
                    return console.error('Error fetching client from pool', err);
                }
                client.query('UPDATE "INFO" SET "' + req.query['column'] + '"=$1 WHERE "sampleHash"=$2', values, function(err, result) {
                    done();
                    if(err) {
                        return console.error('error running query', err);
                    }
                    res.json('OK!');
                });
            });
        } else {
            db.run("UPDATE INFO SET " + req.query['column'] + "=? WHERE sampleHash=?", values, SQLError);
            res.json('OK!');
        }
    } else { res.json('fail.'); }
});

restapi.post('/getData', function(req, res){
    var getDataStart = new Date().getTime();
    parseJSON(req,function(postData) {
        selectWhereGroup = sanitizeQuery(postData)
        if (selectWhereGroup == 'fail') { res.json('fail') }
        else {
            doQuery(postData, selectWhereGroup, res, getDataStart);
        }
    });
});

restapi.get('*', function(req, res) {
    options['path'] = SeQCPath + req.url
    http.get(options, function(getRes) {
        var body = '';
        getRes.on('data', function(chunk) {
            body += chunk;
        });
        getRes.on('end', function() {
            body = body.split('{{hostname}}').join(hostname);               // Replaces all occurences of '{{hostname}}' in AC.GT's code with your servers hostname.
            body = body.split('{{port}}').join(port);                       // Replaces all occurences of '{{port}}' in AC.GT's code with your server's port.
            getRes.headers['content-length'] = body.length;                 // Because reasons.
            res.writeHead(200,getRes.headers);                              // Set the new headers so we dont serve css/etc as html.
            res.write(body);
            res.end();
        });
    }).on('error', function(e) {
        console.log("Got error: " + e.message);
    }); 
});

restapi.listen(port);

console.log("... webservice is up and running. Still checking to make sure we have a connection to your database...");

//'''
