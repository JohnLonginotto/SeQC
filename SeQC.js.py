#!/usr/bin/env python
version = 1//1 ; ''' These two lines are
/*                   really important!  '''

# The first half of this program is written in python. It will count reads in BAM/SAM files using different metrics (per chromosome, per read length, etc).
# To learn more visit http://ac.gt/seqc or run "python ./SeQC.js.py --help"
# The second half of this program is written in JavaScript for use in Node.js. It will create a webserver that can interact with the database created by 
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
import re
import sys
import csv
import json
import time
import types
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
try: import hts ; from cffi import FFI as ffi; htspythonInstalled = True
except: htspythonInstalled = False
try: import psycopg2 ; postgresInstalled = True
except: postgresInstalled = False  

###########################
## Define built-in stats ##
##########################################################################################################
##                                                                                                      ##
## The first thing we do is define our built-in stats. The built-ins are all simple methods to get data ##
## out of the BAM/SAM that already exists. For example chromosome names, template length, CIGAR, etc.   ##
## The names used for each of the built-in stats is exactly the same name as defined by the SAM spec.   ##
## The big exception to all this is FLAG. I refuse to be part of the madness that is talking about SAM  ##
## flags by their hexidecimal sum. Who on earth thought that was a good idea!? I get that is how it is  ##
## stored in the BAM file, but it is akin to talking about DNA as the sum of its 2-bit binary. I'm sure ##
## that even the authors of samtools couldn't tell you what '0x4d3' or its integer equivilent '1235'    ##
## translates to. However, once you know alphabetic notation for flags, 'ABcdEfGHijKl' can be parsed by ##
## humans. Take 2 minutes to look at http://ac.gt/flags and be done with -f and -F forever.             ##
##                                                                                                      ##
## Don't forget - the beauty of SeQC is in how you can easily combine built-in stats and your own code  ##
## to generate new stats in just 2 or 3 lines of python. Reading the SAM/BAM files, parallelization,    ##
## dependency loading, database creation, and vizulization is all taken care of by SeQC. So if you are  ##
## writing new stats, please visit http://ac.gt/seqc/stats for info on how to write stats efficiently.  ##
## It explains how to use dependencies, the "hidden built-ins", and many example stats.                 ##
##                                                                                                      ##
##########################################################################################################
pass

#######################
## Fetch local stats ##
##########################################################################################################
##                                                                                                      ##
## Data that is not already in the BAM/SAM file, such as GC% per read (which has to be calculated, and  ##
## although it might seem obvious, there are actually multiple non-compatible ways to calculate it) are ##
## kept in separate ".stat" files in the same directory as SeQC for simplicity. This is done so people  ##
## can very quickly see exactly how a stat is calculated. External modules take the same format as the  ##
## built-in modules, and it is highly recommended that you use dependencies rather than use pysam or    ##
## hts-python primitives (e.g. don't use line.seq to get the DNA, rather set SEQ as a dependency and it ##
## will be much faster!                                                                                 ##
##                                                                                                      ##
##########################################################################################################

INFO = { 'fileReader': None } # This holds information that a stat module might need to initialize
availableStats = {}           # This is where all the modular stats end up.
def addStat(stat_name,compatible):
    # Do some very basic checks
    if type(stat_name)  is not str:       print '\nERROR: The name of the module from file ' + str(thisFile) + ' is not a string? (this should be the first value for addStat)'; exit()
    if type(stat) is not types.ClassType: print '\nERROR: The module ' + stat_name + ' is not a class? (this should be the second value for addStat)'; exit()
    if type(compatible) is not list:      print '\nERROR: The module ' + stat_name + ' does not provide a lost of other compatible modules. (this should be the third value for addStat)'; exit()
    for md5 in compatible:
        if not re.match(r"[a-f\d]{32}", md5.lower()): print '\nERROR: The MD5 string "' + md5 + '" in module ' + stat_name + ' is not valid MD5.'; exit()
    # Hash the module and init a class object
    with open(thisFile,'rb') as f:
        method = re.sub("addStat.*",'', f.read(), flags=re.MULTILINE|re.DOTALL)
        md5 = hashlib.md5(method)
        if stat_name not in availableStats:
            availableStats[stat_name] = {
                'compatible': json.dumps(compatible),
                'class':      stat, 
                'hash':       md5.hexdigest(),
                'init':       stat(INFO),
                'method':     method,
                'path':       thisFile
            }
        else:
            print '\nERROR: There appears to be two stats with the same name ("' + stat_name + '")\n'
            print availableStats[stat_name]['path'],':'
            print availableStats[stat_name]['method']
            print thisFile,':'
            print method, '\n Choose your favorite one (set the older as compatible if it is), or rename one. This will not effect the MD5.'
            exit()

## Import stat modules from SeQC's directory:
SeQC = os.path.abspath(__file__)                          ## Path to SeQC.js.py itself. Used for spawning new processes.
for potentialStat in os.listdir(os.path.dirname(SeQC)):
        thisFile = os.path.join(os.path.dirname(SeQC), potentialStat)
        if thisFile.endswith('.stat'):
            execfile(thisFile)                            ## I should upgrade this to a proper plugin/module manager, but this works fine for now. 

##########################
## Fetch external stats ##
##########################################################################################################
##                                                                                                      ##
## Get data from the database we're writing to. Perhaps in the future some way to get stats from one    ##
##  database but write data to somewhere else? Like a stat directory db.                                ##
## The danger is that stats can run any code they like (good thing) which could be rm -rf / (bad thing) ##
## and just running SeQC with a stat loaded is enough. The security implications are not so severe if   ##
## stat files are coming from the local directory, but from the internet it could be ...                ##
##                                                                                                      ##
##########################################################################################################
pass

########################
## Check stat modules ##
##########################################################################################################
##                                                                                                      ##
## Check to make sure all the modules are in the right format, have all the required parameters, MD5'd  ##
## etc, etc.                                                                                            ##
##                                                                                                      ##
##########################################################################################################

for statName,statData in sorted(availableStats.items()):
    ## Check stats for required values:
    stat = statData['init']
    for param in [ ('DESCRIPTION',list), ('LINKABLE',bool), ('SQL',str), ('METHOD',str) ]:
    # Uncomment below to try stats without exec()
    #for param in [ ('DESCRIPTION',list), ('LINKABLE',bool), ('SQL',str) ]:
        if not hasattr(stat,param[0]): print '\nERROR: The module ' + statName + ' doesnt have the required parameter ' + param ; exit()
        requiredValue = getattr(stat,param[0])
        if type(requiredValue) is not param[1]:
            if param[0] == 'METHOD' and INFO['fileReader'] == None and requiredValue is None: pass # The only exception.
            else: print '\nERROR: The module ' + statName + ' has parameter ' + param[0] + ' but its not ' + str(param[1]) ; exit()

    # Stat name:
    if len(statName) > 10: print '\nERROR: The name of stat ' + statName + ' is too long! Please keep it to 10 or less characters! :('; exit()

    # DESCRIPTION:
    if len(stat.DESCRIPTION) != 2: print '\nERROR: The DESCRIPTION for stat ' + statName + ' needs two values, an explanation and an output example.'; exit()
    if 60 < len(stat.DESCRIPTION[0]): print '\nERROR: The explination for stat ' + statName + ' needs to be under 70 characters. If this is not possible, put in a URL!'; exit()
    if 30 < len(stat.DESCRIPTION[1]): print '\nERROR: The example output for stat ' + statName + ' needs to be under 30 characters. Put N/A if its way too long, or use ellipses..'; exit()
    if any('%' in s for s in stat.DESCRIPTION): print '\nERROR: The module used to print the --help message cant handle "%" :(\nPlease remove them from your ' + statName + ' module.'; exit()

    # SQL:
    if stat.SQL not in ['TEXT','INT','REAL'] and stat.LINKABLE == True:
        if stat.SQL == 'JSON': print '\nERROR: The stat ' + statName + ' uses JSON as the SQL type, but also claims to be LINKABLE (which it cannot be).'; exit()
        print '\nERROR: The stat ' + statName + ' claims to be a linkable stat, however its .SQL value is not "TEXT", "INT", or "REAL".'; exit()

    # METHOD:
    # If this becomes a function other than an object, this might change by giving it an input value and running it and see if the output is hashable for LINKABLE stated

    # index:
    if hasattr(stat,'index') and stat.index is not False and stat.index is not None:
        if stat.LINKABLE:
            print '\nERROR: The module ' + statName + ' is LINKABLE but also has a value for the index. Indexing for linked stats is done by SeQC. Plese delete the index parameter, or set LINKABLE to False.'; exit()
        else:
            if type(stat.index) is not str: print '\nERROR: The value for self.index in module ' + statName + ' must be a string!'; exit()
            # We could do checks here for things like "CREATE INDEX... ", but actually it might be better to let things be more open for the users.

    # viz:
    if hasattr(stat,'viz') and stat.viz is not False and stat.viz is not None:
        if type(stat.viz) is not list or not (0 < len(stat.viz) < 3): print '\nERROR: The value for viz in module ' + statName + ' must be a list, with 1 or 2 values!'; exit()
        if type(stat.viz[0]) is not str: print '\nERROR: The first value for the viz parameter in module ' + statName + ' is not a string, therefore it is not a name!'; exit()
        if len(stat.viz) is 2 and (type(stat.viz[1]) is not str and type(stat.viz[1]) is not False and stat.viz[1] is not None): 
            print '\nERROR: The code for the viz parameter in module ' + statName + ' must be a string!'; exit()

    # dependencies:
    if hasattr(stat,'dependencies') and stat.dependencies is not None and stat.dependencies is not False:
        if type(stat.dependencies) is not list: print '\nERROR: The value for dependencies in module ' + statName + ' must be None or a list!'; exit()
        for dependency in stat.dependencies:
            if dependency not in availableStats: print '\nERROR: The dependency ' + str(dependency) + ' in module ' + statName + ' is not know to SeQC.'

    # before
    if hasattr(stat,'before') and getattr(stat,'before',None) is not None and stat.before is not False:
        if type(stat.before) is not str: print '\nERROR: The value for self.before in module ' + statName + ' must be None or a string!'; exit() 

    # after
    if hasattr(stat,'after') and getattr(stat,'after') is not None and stat.after is not False:
        if type(stat.after) is not str: print '\nERROR: The value for self.after in module ' + statName + ' must be None or a string!'; exit() 

## Create a dependency graph:
dependencyGraph = {}
class CyclicDependencies(Exception): pass
def sort_dependencies(dependencyGraph):
    post_order = []
    tree_edges = {}
    for fromNode,toNodes in dependencyGraph.items():
        if fromNode not in tree_edges:
            tree_edges[fromNode] = 'root'
            for toNode in toNodes:
                if toNode not in tree_edges:
                    try: post_order += get_posts(fromNode,toNode,tree_edges)
                    except CyclicDependencies as e: print e; exit()
            post_order.append(fromNode)
    return post_order
def get_posts(fromNode,toNode,tree_edges):
    post_order = []
    tree_edges[toNode] = fromNode
    for dependency in dependencyGraph[toNode]:
        if dependency not in tree_edges:
            post_order += get_posts(toNode,dependency,tree_edges)
        else:
            parent = tree_edges[toNode]
            while parent != 'root':
                if parent == dependency:
                    raise CyclicDependencies('\nERROR: Modules ' + dependency + ' and ' + toNode + ' have cyclic dependencies!')
                parent = tree_edges[parent]
    return post_order + [toNode]

for statName in availableStats:
    stat = availableStats[statName]['init']
    if hasattr(stat,'dependencies') and type(stat.dependencies) is list :
        dependencyGraph[statName] = stat.dependencies
    else:
        dependencyGraph[statName] = []
topological_order = sort_dependencies(dependencyGraph) # This is the order dependencies must be run in


helpText = ''
for statName,statData in sorted(availableStats.items()):
    stat = statData['init']
    helpText += statName.ljust(10) + ' | ' + stat.DESCRIPTION[0].ljust(60) + ' | ' + stat.DESCRIPTION[1] + '\n'

## Parse user-supplied command line options:
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
    description="Put in multiple BAM/SAM files, get out an SQL database of read statistics.")
parser.add_argument("--analysis", nargs='+', metavar='', action='append',
    help='''Optional. One or more statistic to gather.
Default is to do "--analysis RNAME TYPE FLAG RGID GC --analysis TLEN" \n
Currently avalible stat names:
''' + helpText + '\n')
parser.add_argument("--input", nargs='+', metavar='',
    help='Required. One or more SAM (and/or BAM) files to analyse.')
parser.add_argument("--output", default='myProject.SeQC', metavar='',
    help='Optional. Name of output database (or path for SQLite). Default is "myProject".')
parser.add_argument("--samtools", default='samtools', metavar='', 
    help="Optional, unless the path to samtools is needed and cannot be found.")
parser.add_argument("--quiet", action='store_true',
    help='Optional. No status bars in output. Good for logs, bad for humans.')
parser.add_argument('--writeover', action='store_true',
    help="Optional. Will write over old data even if inputs/analyses are identical")
parser.add_argument("--pguser", metavar='',
    help="Optional, unless you want to use Postgres instead of SQLite. User account name.",)
parser.add_argument("--pgpass", metavar='',
    help="Optional. Password of --pguser.")
parser.add_argument("--pghost", default="localhost", metavar='',
    help="Optional. Hostname of postgres database. Default is localhost.")
parser.add_argument("--cpu", default=2, metavar='', type=int,
    help="Optional. Number of processes/cores you want to use. Default is 2.")
parser.add_argument("--md5", action='store_true',
    help="Optional. Returns the MD5 checksums of all loaded modules (as SeQC sees them) and exits.")
parser.add_argument("--debug", action='store_true',
    help="Required. Probably.")
parser.add_argument("--SAM", help=argparse.SUPPRESS)   # Used internally to tell subprocesses we are reading SAM. Set to either 'stdin' or 'file'.
parser.add_argument("--BAM", help=argparse.SUPPRESS)   # Used internally to tell subprocesses we reading directly via pysam or htspython.
args = parser.parse_args()

## Print MD5s (without anything after the addStat) to make registering other modules as compatible easier.
if args.md5:
    for stat_name,stat in availableStats.items():
        print 'Stat Name:'.rjust(16), stat_name
        print 'Hash:'.rjust(16),stat['hash']
        print 'Path:'.rjust(16),stat['path']
        print 'Compatible With:'.rjust(16),stat['compatible'],'\n';
    exit()

## It takes your stats and returns all stats (so, includes the dependancies, and the dependancies of dependancies..)
def get_all_modules_to_run(stat_name):
    if hasattr(availableStats[stat_name]['init'],'dependencies'):
        immediate_dependencies = set(availableStats[stat_name]['init'].dependencies)
        more_dependencies = set()
        for dependency in immediate_dependencies:
            more_dependencies.update(get_all_modules_to_run(dependency))
        immediate_dependencies.add(stat_name)
        return tuple(immediate_dependencies | more_dependencies)
    else: return (stat_name,)

## "Sanitize" user analysis arguments:
if args.analysis is None:
    args.analysis = [('RNAME', 'TYPE', 'FLAG', 'GC'),('TLEN',)] # Default analyses

if not args.SAM and not args.BAM: 
    print '   [ ' + str(len(availableStats)) + ' Modules Loaded! ]'
# We make each linked group of stats a set, before converting to a tuple, in the event that the user adds the same stat twice, eg. -a gc gc becomes just (gc).
# We also sort by stat name in each linked group of stats so "-a tlen gc" becomes (gc,tlen).
# We then put these tuples into a set, so that "-a gc tlen -a tlen gc" would become just ((gc,tlen)), as a form of redundancy checking.
# Finally we then put all these sorted tuples into a sorted list. ~ phew ~
allGroups = set()
allAnalyses = set()
for group in args.analysis:
    group = set(group)
    for stat in group:
        if stat not in availableStats.keys():
            print '\nERROR: I do not know how to calculate the statistic: ' + stat
            print 'If the .stat file for this statistic is not in the same directory as SeQC, you may need to download it from http://ac.gt/seqc\n'; exit()
        if availableStats[stat]['init'].LINKABLE == False and len(group) != 1:
            print '\nERROR: Unlinkable stats must go in their own --analysis group! Please take "' + stat + '" out of --analysis ' + ' '.join(group); exit()
        allAnalyses.update(get_all_modules_to_run(stat)) # Recurses all dependancies too.
    allGroups.add(tuple(sorted(group)))                  # sorted returns a list but we need a hashable tuple.
args.analysis = sorted(allGroups)                        # A sorted list of sorted tuples.
sorted_analyses = []                                     # This is a non-redundant list of the analyses used (and there dependancies), in the order they need to be run. 
for analysis in topological_order:
    if analysis in allAnalyses:
        sorted_analyses.append(analysis)

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
If this is your first time using SeQC, try "SeQC.js.py --help" or visit http://ac.gt/seqc/ for usage infomation.
'''; exit()

    ## Determine if the user is writing to a SQLite or Postgres database, and determine that the provided details work:
    if args.pguser == None:
        print '   [ Using SQLite for output ]'
        if os.path.isdir(args.output): print '''
ERROR: You have provided an existing DIRECTORY for your output, but I need a file name!
Please use the exact name of the output file you want (you can always rename it later)'''; exit()

        if os.path.isfile(args.output): print '   [ Output file already exists ]'
        else:                           print '   [ Output file does not exist ]'
    else:
        print '   [ Using Postgres for output ]'
        if postgresInstalled:
            if args.pgpass == None:
                if args.debug: 
                    password = 'PASSWORD_HIDDEN' # we dont want --debug to print passwords! Since --debug on the parent never executes jobs anyway, this is fine.
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
            print '\nERROR: The input file ' + str(inputFile) + ' could not be accessed. Are you sure it exists and we have permissions to read it? (continuing without it)'
        else: usedInputs.append(inputFile)

    ## Most parameters to all subprocesses are static and never change, such as the database connection details or --writeover. 
    ## Here we collect them all into 1 string called explicitStats:
    explicitStats = ''
    for group in args.analysis: explicitStats += ' --analysis ' + ' '.join(group)
    if args.writeover: explicitStats += ' --writeover'
    if args.pguser != None:
        explicitStats += ' --pguser "' + args.pguser + '" --pgpass "' + password + '" --pghost "' + args.pghost + '"'

    ## This function checks to see if a file is binary (BAM) or ASCII (SAM) by looking for a "null byte" (0x00) in the first 20Mb of the file.
    def bamCheck(fileName):
        with open(fileName, 'rb') as xamfile:
            for i in range(10): # 10 tries to find a null byte in 2Mb chunks.
                section = xamfile.read(2048)
                if '\0' in section:            return True  # BAM file.
                if len(section) < 2048:
                    if i == 0 and not section: return None  # Empty file.
                    else:                      return False # SAM file.

    ## We only need samtools if we have one or more BAM files. Here we check all files with the bamCheck function, 
    ## and if a binary file is found, we check if we can find/execute samtools.
    gotSAM,gotBAM,samtoolsInstalled = False, False, None # If user gives lots of files, we could have both.
    DEVNULL = open(os.devnull, 'wb')
    if args.debug: STDERR = subprocess.STDOUT
    else: STDERR = DEVNULL
    for inFile in usedInputs:
        if bamCheck(inFile):
            gotBAM = True
            if not args.samtools:
                exitcode = subprocess.call('samtools', stdout=DEVNULL, stderr=DEVNULL, shell=True)
                if exitcode != 1: samtoolsInstalled = False # Samtools is weird, in that calling it with no parameters returns exitcode 1 rather than 0.
                else: args.samtools = 'samtools'; samtoolsInstalled = True
            else:
                code = subprocess.call(args.samtools, stdout=DEVNULL, stderr=DEVNULL, shell=True)
                if code != 1:
                    print '\nERROR: You have tried to calculate statistics on a BAM file, but the path to samtools you have provided:'
                    print args.samtools,' does not work :('
                    exit()
                else: samtoolsInstalled = True

    ## Check that all the used stats have a method for this kind of data:
    blocking = {'sam':False,'pysam':False,'htspython':False}
    for stat in sorted_analyses:
        if gotSAM:
            # Here we can only use the sam fileReader. Perhaps I should extend this to allow for pysam/hts if they can work on SAM files... hm..
            availableStats[stat]['init'] =  availableStats[stat]['class']({'fileReader':'sam'}) # re-init with the sam fileReader
            if availableStats[stat]['init'].METHOD == None: blocking['sam'] = stat
        if gotBAM:
            if htspythonInstalled:
                availableStats[stat]['init'] = availableStats[stat]['class']({'fileReader':'htspython'}) # re-init with htspython
                if availableStats[stat]['init'].METHOD == None: blocking['htspython'] = stat # Since we cannot use it.
            else: blocking['htspython'] = stat
            if pysamInstalled:
                availableStats[stat]['init'] = availableStats[stat]['class']({'fileReader':'pysam'}) # re-init with pysam
                if availableStats[stat]['init'].METHOD == None: blocking['pysam'] = stat
            else: blocking['pysam'] = stat
            if samtoolsInstalled:
                availableStats[stat]['init'] = availableStats[stat]['class']({'fileReader':'sam'}) # re-init with sam
                if availableStats[stat]['init'].METHOD == None: blocking['sam'] = stat
            else:blocking['sam'] = stat
    if all(blocking.values()):
        print '\nERROR: The combination of stats you have selected can not be gathered, because there is no method for '
        print 'reading the file (sam/pysam/htspython) that all the stats can use, specifically:'
        for x,y in blocking.items(): print y,'does not work with',x
        print '\nEither hack on the stat modules to get them to work for these other file readers, or if you are not linking these conflicting stats, run SeQC multiple times.\n'
        exit()
    elif not blocking['htspython']: INFO = {'fileReader':'htspython'}
    elif not blocking['pysam']:     INFO = {'fileReader':'pysam'}
    elif not blocking['sam']: INFO = {'fileReader':'sam'}
    print '   [    Using',INFO['fileReader'],'   ]'

    ## The table schema for the INFO table. 
    ## The INFO table stores metadata on all the analysed files, as well as which stats (or linked groups of stats) have been collected for each file.
    ## It is basically 1 row per analysed BAM/SAM file (or rather, per input file MD5 checksum).
    ## Regarding data privacy, SeQC make no distinction between users accessing data - be it in-house, or across the globe. It also makes no distinction between
    ## users connecting directly, and users finding your database through the public listings (if you decide to make your data publicly available).
    ## SeQC does allow you to hide any INFO column from users via SERVER settings. file name and path for example are quite popular ones to hide. You cannot, however,
    ## hide *some* of the data in json_stats. It's all or nothing. This is why any data going into json_stats is always gathered explicitly.
    infoTableCreation = '''CREATE TABLE "INFO" (
                            "hash"          TEXT PRIMARY KEY,
                            "path"          TEXT,
                            "size"          TEXT,
                            "header"        TEXT,
                            "analyses"      TEXT,
                            "file_name"     TEXT,
                            "json_stats"    TEXT,
                            "last_updated"  TEXT,
                            "creation_time" TEXT,
                            "total_reads"   INT
                        )'''
    ## The data gets added to INFO by subprocesses (after all stats for that input file have been calculated).

    ## The table schema for the SERVER table. 
    ## Very simple key-value storage. Will contain information on what the landing page should look like, who can access the webserver, etc.
    ## Data here should never leave the server itself. If it already exists in the output database we do nothing.
    serverTableCreation = '''CREATE TABLE "SERVER" (
                            "key"   TEXT PRIMARY KEY,
                            "value" TEXT
                        )'''
    serverTableData   = [
                            ('public'   , 'False'    ),
                            ('listen_on', '8080'     ),
                            ('bind_to'  , '127.0.0.1')
                        ]

    ## The table schema for the METHOD table. It stores the code for both stats AND visualization code.
    ## Regarding stat data, in order to make sure different users can write their own stats and share/compare the results with one another, we need a
    ## way to make sure that the stats calculated differently are not compared side-by-side without some sort of warning (even if they have the same name). 
    ## We do this by MD5ing the code that generates the stat (but after cutting out it's self.COMPATIBLE data and sticking it in the compatible column).
    ## We do essentially the same thing with the visualization code, except this code does not have compatibility issues.
    ## Please bear in mind that if you share your data, you also share the methods to create it. 
    methodTableCreation = '''CREATE TABLE "METHOD" (
                            "hash"       TEXT PRIMARY KEY,
                            "type"       TEXT,
                            "name"       TEXT,
                            "code"       TEXT,
                            "compatible" TEXT
                        )'''
    methodTableData = []
    for statName,stat in availableStats.items():
        if statName in sorted_analyses:
            methodTableData.append((stat['hash'],'stat',statName,stat['method'],stat['compatible']))

    ''' JOHN: you might be able to merge sqlite and postgres code here '''

    ## Create the INFO table (if not already present) for SQLite:
    if args.pguser == None: # No pguser means SQLite!
        try:
            con = sqlite3.connect(args.output, timeout=120)
            cur = con.cursor()
            cur.execute("SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='INFO';")
            if cur.fetchone()[0] != 1:
                print '   [    Creating INFO table    ]'
                cur.execute(infoTableCreation)
            cur.execute("SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='SERVER';")
            if cur.fetchone()[0] != 1:
                print '   [   Creating SERVER table   ]'
                cur.execute(serverTableCreation)
                cur.executemany('INSERT INTO "SERVER"(key,value) values (?,?)',serverTableData)
            else:
                pass # Currently, if there is already a SERVER table, we just don't touch it. Users can do what they like.
            cur.execute("SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='METHOD';")
            if cur.fetchone()[0] != 1:
                print '   [   Creating METHOD table   ]'
                cur.execute(methodTableCreation)
                cur.executemany('INSERT INTO "METHOD"(hash,type,name,code,compatible) values (?,?,?,?,?)', methodTableData)
                print '   [     All tables created    ]'
            else:
                print '   [   Updating METHOD table   ]'
                for methodHash,methodType,methodName,methodCode,methodCompatible in methodTableData:
                    cur.execute('SELECT compatible FROM "METHOD" WHERE hash=?',(methodHash,))
                    result = cur.fetchone()
                    if result:
                        methodCompatible = json.dumps(list(set(json.loads(methodCompatible) + json.loads(result[0])))) # merge compatibility lists
                    cur.execute('INSERT OR REPLACE INTO "METHOD"(hash,type,name,code,compatible) values (?,?,?,?,?)',(methodHash,methodType,methodName,methodCode,methodCompatible))
            con.commit() # So we can roll all of the above back if there's an error.
        except sqlite3.Error, e:
            print '''
ERROR: Something went wrong creating the output database. The exact error was:
''' + str(e) + '''
If you do not know what caused this error, please e-mail john@john.uk.com with the error message and I will help you :)'''; exit()
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
                print '   [ Creating Postgre database ]'
                cur.execute('CREATE DATABASE "' + args.output + '"')
                con.commit()
            cur.close();con.close()
            con = psycopg2.connect(dbname=args.output,user=args.pguser,host=args.hostname,password=password) # ...and connect to it.
            cur = con.cursor()
            cur.execute("SELECT * from information_schema.tables where table_name='INFO'")
            if cur.rowcount != 1:
                print '   [    Creating INFO table    ]'
                cur.execute(infoTableCreation)
            cur.execute("SELECT * from information_schema.tables where table_name='SERVER'")
            if cur.rowcount != 1:
                print '   [   Creating SERVER table   ]'
                cur.execute(serverTableCreation)
                cur.executemany('INSERT INTO "SERVER"(key,value) values (%s,%s)',serverTableData)
            cur.execute("SELECT * from information_schema.tables where table_name='METHOD'")
            if cur.rowcount != 1:
                print '   [   Creating METHOD table   ]'
                cur.execute(methodTableCreation)
                cur.executemany('INSERT INTO "METHOD"(hash,type,name,code,compatible) values (?,?,?,?,?)', methodTableData)
                print '   [     All tables created    ]'
            else:
                print '   [   Updating METHOD table   ]'
                for methodHash,methodType,methodName,methodCode,methodCompatible in methodTableData:
                    cur.execute('SELECT compatible FROM "METHOD" WHERE hash=%s',(methodHash,))
                    result = cur.fetchone()
                    if result:
                        methodCompatible = json.dumps(list(set(json.loads(methodCompatible) + json.loads(result)))) # merge compatibility lists. Also postgres returns a string not tuple so no result[0]
                        cur.execute('INSERT INTO "METHOD"(hash,type,name,code,compatible) VALUES(%s,%s,%s,%s,%s) ON CONFLICT (hash) DO UPDATE SET compatible = EXCLUDED.compatible',
                            (methodHash,methodType,methodName,methodCode,methodCompatible))
            con.commit() # But apparently I'M the one who can't commit..
            cur.close();con.close()
        except psycopg2.ProgrammingError, e:
            print '''
ERROR: Something went wrong creating the output database. The exact error was:
''' + str(e) + '''
If you do not know what caused this error, please e-mail longinotto@immunbio.mpg.de with the error message and I will help you :)'''; exit()
        finally:
            if cur: cur.close()
            if con: con.close()

    ## This function fires off the subprocesses which actually analyse the input BAM/SAM data.
    ## It will be called in a loop later as we process the input files.
    def doSeQC(inputFile):
        if bamCheck(inputFile):
            if INFO['fileReader'] == 'htspython' or INFO['fileReader'] == 'pysam':
                subprocessCommand = 'python "' + SeQC + '" --input "' + inputFile + '" --BAM ' + INFO['fileReader'] + ' --output "' + args.output + '" ' + explicitStats
                if args.debug: print subprocessCommand # Don't run, just print what would have run.
                else: return subprocess.Popen(subprocessCommand, stdout=subprocess.PIPE, stderr=STDERR, shell=True, executable='/bin/bash')
            else:
                subprocessCommand = '"' + args.samtools + '" view "' + inputFile + '" | python "' +  SeQC + '" --input "' + inputFile + '" --output "' + args.output + '" --SAM stdin ' + explicitStats
                if args.debug: print subprocessCommand
                else: return subprocess.Popen(subprocessCommand, stdout=subprocess.PIPE, stderr=STDERR, shell=True, executable='/bin/bash')
        else:
            subprocessCommand = 'python "' + SeQC + '" --input "' + inputFile + '" --SAM file --output "' + args.output + '" ' + explicitStats
            if args.debug: print subprocessCommand
            else: return subprocess.Popen(subprocessCommand, stdout=subprocess.PIPE, stderr=STDERR, shell=True, executable='/bin/bash')

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
            except: tidy('\nERROR: Process analysing ' + theFile + ' unexpectedly stopped?!')

            if out in ['0','1','2','3','4','5','6','7','8','9']: outputs[theFile][1] += int(out)
            elif out == '#': outputs[theFile] = ['Calculating MD5         |',0]                                     ## Statuses are always
            elif out == '|': outputs[theFile] = ['Calculating Statistics  |',0]                                     ## exactly 25 characters
            elif out == '&': outputs[theFile] = ['Converting to SQL table |',0]                                     ## in width and include
            elif out == '@': outputs[theFile] = ['Writing To Database     |',0]                                     ## the starting pipe "|"
            elif out == '$': outputs[theFile] = ['Postgres Table Indexing |',0]
            elif out == '%': tidy('File ' + theFile + ' skipped as it is already in the database.')
            elif out == '?': tidy('DATA ERROR: ' + theFile + ' aborted by SeQC. Rerun with --debug for more info.')
            elif out == '!': tidy('Completed: ' + theFile)
            else:
                print 'x'+out+'x'
                tidy('SeQC ERROR: ' + theFile + ' failed, likely due to a bug in SeQC. Run with "--debug" for more info.')

        ## Print some beautiful status bars ;)
        progressBarSpace = int(terminalCols/2)                              ## Use half the screen for the progress bar and status.
        fileNameSpace = terminalCols - progressBarSpace                     ## Use the rest for the file name (full path).
        progressScaling = (progressBarSpace -27) / 100.                     ## -27 because we dont want the scaling factor to know about the status, arrow, or final pipe.
        arrow = '' if arrow else '>'                                        ## '' is falsey, while '>' is truthy, so this toggles it.
        for theFile, storedOut in outputs.items():
            if len(theFile) > fileNameSpace:
                left = '<-- ' + (theFile + ': ')[-fileNameSpace+4:]         ## <-- to indicate the file name was truncated
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
    total_reads, file_hash = MD5andCount(inputFile,args.samtools)

    ## Now we check if any analyses have been performed on this sample before.
    ## If they have and --writeover is not set, we dont do the analysis again.
    if args.pguser == None: con = sqlite3.connect(args.output, timeout=120)
    else:                   con = psycopg2.connect(dbname=args.output, user=args.pguser, host=args.pghost, password=args.pgpass)
    cur = con.cursor()
    ## We have two options to find out if the table exists - check the SQL database itself, or check the INFO table. Which is 'more correct'? I think checking the INFO table
    ## is, as if it's not there, it doesn't matter whats already in the database, we're expected to add it in over the top - but below is the code to check the database directly
    ## cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='" + file_hash + "';")                                                            ## SQLite
    ## cur.execute("SELECT 1 FROM pg_catalog.pg_class WHERE relkind = 'r' AND relname = '" + args.output + "' AND pg_catalog.pg_table_is_visible(oid) LIMIT 1") ## Postgres
    cur.execute("SELECT \"analyses\", \"json_stats\" FROM INFO WHERE \"hash\"='" + file_hash + "';")
    result = cur.fetchone()
    # "analyses" is a dictionary (or a JSON object...) where the key is the --analysis group after a '_'.join(group), with all sorts of meta data as values.
    # Because analysis groups are always stored/requested in alphabetical order, the key for any given group of analyses is easy to work out.
    if result != None:
        existing_analyses = json.loads(result[0])
        existing_json_stats = json.loads(result[1])
        if type(existing_analyses) != dict:
            if args.debug: print "\nERROR: In your output database's INFO table, the analyses column for the file " + args.input + ' (' + file_hash + ') is not a dictionary/object?!'; exit()
            else: print '?'
            exit()
        if not args.writeover:
            skip = set()
            for group in args.analysis:
                if '_'.join(group) in existing_analyses:
                    skip.add(group)
            for x in skip:
                args.analysis.remove(x)
            if len(args.analysis) == 0:
                sys.stdout.write('%'); sys.stdout.flush() # Skipping
                exit()
    else:
        existing_analyses = {}
        existing_json_stats = {}
    if cur: cur.close()
    if con: con.close()

    ## Finally, we now start analyzing the data...
    ping.change('|',total_reads)
    data = [collections.defaultdict(int) for x in range(0,len(args.analysis))] # a list of dictionaries, one for every analysis group.
    if args.SAM:
        # for stat in availableStats:
        #     availableStats[stat].process = availableStats[stat].SAM
        # if args.SAM == 'stdin':
        #     sys.argv = '' ## Done so fileinput takes stdin and not args.
        #     inputData = csv.reader(fileinput.input(), delimiter='\t')
        # else: 
        #     inputData = csv.reader(fileinput.input(inputFile), delimiter='\t') # fileinput also can read a file on disk
        # header = json.dumps('')
        pass
    elif args.BAM:
        for stat in sorted_analyses:
            availableStats[stat]['init'] = availableStats[stat]['class']({'fileReader':args.BAM})
        if args.BAM == 'htspython':
            inputData = hts.Bam(inputFile, "rb")
            header = json.dumps('htspython currently does not support reading the BAM header. If this upsets you please mail the author Brent Pedersen via GitHub:)')
        elif args.BAM == 'pysam':
            inputData = pysam.Samfile(inputFile, "rb")
            header = json.dumps(inputData.header)

    ## Do .before work - run ONCE per module.
    for analysis in sorted_analyses:
        if hasattr(availableStats[analysis]['init'],'before'):
            exec(availableStats[analysis]['init'].before)

    ## CPython - the program that usually runs your python scripts - sucks at inlining function calls.
    ## Pypy, a sepurate program that can execute python scripts, doesn't have this problem - however from a recent
    ## poll in biostars (https://www.biostars.org/p/181266/) it appears few people have pypy installed.
    ## As a result, i'm going to inline the code at run-time below, which is why stat functions must be strings not actual functions:
    execute  = '\n'
    execute += 'for reads_processed, read in enumerate(inputData):\n'
    execute += '    ping.pong(reads_processed)\n'
    for analysis in sorted_analyses:
        analysis_method = availableStats[analysis]['init'].METHOD.splitlines()
        whitespace_count = 0
        whitespace_marker = None
        for line in analysis_method:
            whitespace_count = len(line) - len(line.lstrip())
            if whitespace_count: # if there is indentation.
                whitespace_marker = line[0]
                break
        if whitespace_marker is None:
            reindended_method = ''
            for line in analysis_method:
                reindended_method += '    ' + line + '\n'
        else:
            reindended_method = ''
            for line in analysis_method:
                indentation_level = 0
                while line[:whitespace_count] == whitespace_marker*whitespace_count:
                    line = line[whitespace_count:]
                    indentation_level += 1
                reindended_method += ((1+indentation_level)*'    ') + line + '\n' # +1 to make sure everything is indented at least once for our "for reads_processed...." loop.
        execute += reindended_method

    for idx, analysis in enumerate(args.analysis):
        if all([ availableStats[x]['init'].LINKABLE for x in analysis ]):
            if len(analysis) == 1:
                execute += '    data[' + str(idx) + '][('+ analysis[0] + ',)] += 1\n'
            else: 
                execute += '    data[' + str(idx) + '][('+ ','.join(analysis) + ')] += 1\n'
    if args.debug: print execute
    exec(execute)

    ## Set up some generic database-specific things:
    if args.pguser == None:
        con = sqlite3.connect(args.output, timeout=120)
        con.isolation_level = 'EXCLUSIVE'
        con.execute('BEGIN EXCLUSIVE')
        cur = con.cursor()
        delim = '?'
    else:
        con = psycopg2.connect(dbname=args.output, user=args.pguser, host=args.pghost, password=args.pgpass)
        cur = con.cursor()
        delim = '%s'
    new_json_stats = {} # stores the json_stats data
    new_analyses = {} # stores the analysis metadata

    ## From here on out, we process things one analysis group at a time, until its time to update the INFO table.
    for group_idx, group in enumerate(args.analysis):
        tableName = file_hash + '_' + '_'.join(group)            # e.g. 2dfdd09ddb42ef2fab575012815a06b6_TLEN_RNAME. Not used for JSON.
        cur.execute('DROP TABLE IF EXISTS "' + tableName + '"') # if --writeover wasnt set, and there was already data, the analysis group would have been dropped earlier.
        cur.execute('DROP INDEX IF EXISTS "' + tableName + '_INDEX"')
        ping.change('&',len(args.analysis))

        ## Now what we do depends on if we're linkable or not:
        if all([ availableStats[x]['init'].LINKABLE for x in group ]):
            ## data is currently a list of dictionaries, where the keys of those dictionaries are immutable tuples and the values are read counts. This is not a bag of fun.
            ## So for every dictionary we convert it to a mutable list-of-lists (i.e. a table, for the SQL database). We call this table "table".
            table = []
            for stat_columns,read_count in data[group_idx].items():
                table.append(list(stat_columns) + [read_count])
                # ping.pong(group_idx) # Could turn this on if this convertion consistently takes longer than 1 second ...

            ## Do .after work. As you can see, this happens for every time the analysis is used in an analysis group, not just once like .before
            ## This is because .after is often used to modify the generated data, like for RNAME, it converts numerical chromosome ids to real names.
            for column, analysis in enumerate(group):
                if hasattr(availableStats[analysis]['init'],'after'):
                    ping.pong(group_idx)
                    exec(availableStats[analysis]['init'].after)

            ## Dont add anything when there's nothing to add:
            if len(table) == 0 or (len(table) == 1 and all(value == None for value in table[0])):
                if debug: print '\nERROR: No data was collected from the (' + ' '.join(group) + ') analysis, so nothing at all was written to the database!'
                continue # Move on to next analysis group

            ## OK, now time to add that data to a new table.
            ping.change('@',len(args.analysis))

            ## Set up LINKABLE=True specific SQL stuff:
            if args.pguser == None:
                createColumn = ', '.join([stat+' '+availableStats[stat]['init'].SQL for stat in group]) + ', counts INTEGER'
                columnQ      = ', '.join([delim]*(len(group)+1)) # +1 for the counts column
            else:
                pgMegaString = '' # For postgres we do something a bit funny - we use it's COPY function and provide a huge string, as it's really fast.
                column_idxs = range(len(table[0]))
                for row_idx in xrange(0,len(table)):
                    for column_idx in column_idxs:
                        if table[row_idx][column_idx] == None:
                            table[row_idx][column_idx] = '\N' ## postgres identifies null as \N  -  note that its not a special symbol, it's literally just a backslash and a capital N. This can be changed with " WITH NULL AS 'tumblr' " if you'd rather 'tumblr' was the signifier for something with no value.
                    pgMegaString = '\t'.join(table[row_idx]) + '\n'

            ## Create table
            cur.execute('CREATE TABLE "' + tableName + '"(' + createColumn + ')')
            if args.pguser == None:
                cur.executemany('INSERT INTO "' + tableName + '" VALUES ( ' + columnQ + ' )', table)
            else:
                cur.copy_from(StringIO.StringIO(pgMegaString), '"' + tableName + '"')
            con.commit()
            new_analyses['_'.join(group)] = {
                'rows'       :   len(table),
                'stat_hashes':   [ availableStats[stat]['hash'] for stat in group ],
                'performed'  :   int(time.mktime(time.gmtime()))*1000 # milliseconds since unix epoch in GMT/UMT
            }
        else:
            analysis = group[0]
            exec('unlinkable_result = '+analysis) # I don't want to force unlinkables to have to put their data anywhere special other than their variable name, and this works.
            if availableStats[analysis]['init'].SQL == 'JSON':
                try: json.dumps(unlinkable_result)
                except TypeError:
                    print '\nERROR: The module "' + analysis + ' has returned ' + str(type(unlinkable_result)) + ' data that can not be JSON encoded.'
                    print '       This is a bug in the module. Please advise whomever wrote this stat that it occured and how (if you can)'; exit()
                new_json_stats[analysis] = unlinkable_result
                new_analyses['_'.join(group)] = {
                    'rows'       :   None,
                    'stat_hashes':   [ availableStats[stat]['hash'] for stat in group ],
                    'performed'  :   int(time.mktime(time.gmtime()))*1000 # milliseconds since unix epoch in GMT/UMT
                }
            else:
                if type(unlinkable_result) != list:
                    print '\nERROR: The module ' + analysis + ' has returned a ' + str(type(unlinkable_result)) + ' and not a list of lists.'
                    print '       This is a bug in the module. Please advise whomever wrote this stat that it occured and how (if you can)'; exit()
                for row in unlinkable_result:
                    if type(row) != list and type(row) != tuple:
                        print '\nERROR: The module ' + analysis + ' has returned a list of ' + str(type(unlinkable_result)) + ' and not a list of lists.'
                        print '       This is a bug in the module. Please advise whomever wrote this stat that it occured and how (if you can)'; exit()

                createColumn = ', '.join([x[0]+' '+x[1] for x in availableStats[analysis]['init'].SQL]) # tuple/list of (column_name,column_type) tuples/lists
                cur.execute('CREATE TABLE "' + tableName + '"(' + createColumn + ')')
                # Here we won't use Postgres' COPY function, since this gives the user a bit more flexibility in adding special postgres-specific datatypes (lists/etc).
                cur.executemany('INSERT INTO "' + tableName + '"(' + ', '.join([x[0] for x in availableStats[analysis]['init'].SQL]) + ') values ( ' + ', '.join([delim for x in availableStats[analysis]['init'].SQL]) + ' )', unlinkable_result)
                con.commit()
                new_analyses['_'.join(group)] = {
                    'rows'       :   len(unlinkable_result),
                    'stat_hashes':   [ availableStats[stat]['hash'] for stat in group ],
                    'performed'  :   int(time.mktime(time.gmtime()))*1000 # milliseconds since unix epoch in GMT/UMT
                }

    ## DONE MAKING TABLES! Now we just have to update the INFO table!
    ## First we grab some meta data about the file
    file_name = os.path.basename(args.input[0])                  # Always replaces
    file_path = os.path.abspath(args.input[0])[:-len(file_name)]  # Always replaces
    creation_time = int(time.mktime(time.gmtime(os.path.getmtime(args.input[0]))))*1000 # NEVER replaces. Creation of file, again milliseconds from unix epoch in GMT/UMT
    analysis_time = int(time.mktime(time.gmtime()))*1000    # Always replaces. Time right now.
    sample_size = os.path.getsize(args.input[0])                 # If this is different, you're going to get an error (same hash different file size? hm...)

    ## If this is the first time this sample has had an entry in the INFO table, we have to make a new row:
    if len(existing_analyses) == 0:
        for key,value in new_analyses.items():   existing_analyses[key] = value
        for key,value in new_json_stats.items(): existing_json_stats[key] = value
        existing_analyses   = json.dumps(existing_analyses)
        existing_json_stats = json.dumps(existing_json_stats)
        cur.execute('INSERT INTO "INFO" ("hash",    "path",    "size",      "header", "analyses",        "file_name", "json_stats",        "last_updated", "creation_time", "total_reads") VALUES (' +  ', '.join((delim,)*10) + ')',
                                        ( file_hash, file_path, sample_size, header,   existing_analyses, file_name,   existing_json_stats, analysis_time,  creation_time,   total_reads )                                        )
    ## An INFO row already exists for this sample...
    else:
        for key,value in new_analyses.items():   existing_analyses[key] = value
        for key,value in new_json_stats.items(): existing_json_stats[key] = value
        existing_analyses   = json.dumps(existing_analyses)
        existing_json_stats = json.dumps(existing_json_stats)
        sql = 'UPDATE "INFO" SET "path"=$, "analyses"=$, "file_name"=$, "json_stats"=$, "last_updated"=$ WHERE "hash"=$'.replace('$',delim)
        cur.execute(sql, (file_path, existing_analyses, file_name, existing_json_stats, analysis_time, file_hash))
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
