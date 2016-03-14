# Hello :) 
# I will explain the formatting for modules by example :) First some do's and dont's:

# A lot of the code is boilerplate. Copy and paste another module to get started :)
# Anything in CAPS is required, and anything in lowercase is optional.
# It is always better to use another module as a dependency than fetch the value yourself, because you may otherwise end up fetching values twice.
# Try not to use code that is sam/pysam/htspython specific if you can help it.
# The code in this file is executed in-place, so it has full access to the namespace/data it runs in. In short, that is a lot of power, which is good, but try 
# not to write-over existing variable names or your modules won't work. I like to prefix my variables with the stat name first, like "FLAG_intToFlag".
# self.METHOD, self.before and self.after have to be strings which will get exec()uted. This is not ideal, but it's significantly faster than the alternative.
# Once you write a stat, it may take on a life of its own. It will be shared with the data that it creates if you choose to make that data public. 
# For this reason, consider popping some contact details into a comment at the top to let other people know who put in the hard work ;-)

class stat:                   # dont change this
    def __init__(self,INFO):  # or this.
        self.DESCRIPTION =  ['The ID of the read', 'HWI-ST51_039:8:2:152:511/1'] # Description and example result. Keep both short or you'll get an error.

        self.LINKABLE   = True                  # What does it mean to be linkable?
                                                # SeQC was designed to answer questions that start with:
                                                # "How many reads have/are/contain ..." or "What is the distribution of  ... for all my reads.", for example:
                                                # "How many reads are on Chromosome 1"  or "What is the distribution of TLEN for all my reads."
                                                # For this reason, all stats that SeQC produces that are linkable have a column called "reads" and the
                                                # sum of that column adds up to the total number of reads in your sample file. This way many stats can be
                                                # linked together. E.g. "chromosomes and duplicated linked together" might look like:
                                                # rname      duplicated      reads
                                                # chr1       True            234334           [note that there is no stat called "duplicated", but you could make one easily]
                                                # chr1       False           3298272
                                                # chr2       True            ...
                                                # So if a stat is linkable in this way, self.LINKABLE must be True.

                                                # Stats that are not linkable, are stats which do not follow the above description. For example:
                                                # "The total number of mapped base pairs" or "The distribution of sequencing depth".
                                                # Often unlinkable stats can only be calculated after all the data has been read, like in the first example.
                                                # Some unlinked stats may not even be producable by SeQC, but can still be stored in the database and shared/visualized.
                                                # Unlinked stats broadly fall into two catagories:
                                                # 1) stats which do not benefit from being stored in SQL, such as single values or non-tabluar data.
                                                # 2) stats which do benefit from SQL and indexing.
                                                # The former will be stored in JSON in the 'json_stats' column of the sample's row in the INFO table.
                                                # The latter will be made into its own table (like a linked data table), but you must also specify the schema, indexes, and viz.

        self.SQL = 'TEXT'                       # Tell SeQC the datatype of the stat when stored in the SQL database.
                                                # For "self.LINKABLE == True" data, there are three options:
                                                # For data that is a string (and therefore catagorical data) use 'TEXT',
                                                # e.g. for "how many reads are on Chromosome 1" we use the RNAME stat, which stores chromosome names and uses 'TEXT' 
                                                # For data that is an integer or decimal float, use "INT" or "REAL" respectively (these are SQL terms by SQLite and Postgres).
                                                # e.g. "what is the distribution of TLEN?" or "What is the distribution of GC%?"

                                                # For data that is "self.LINKABLE == False", there are two options:
                                                # 1) Use "JSON" if you wish to store the stat in the INFO table with the sample metadata, in the "json_stats" column.
                                                # It will be stored in a JSON object where the key is the stat name, and the value is your str/int/list/tuple/dict data.
                                                # If the value is not JSON, it must be a list of (column name/column type) tuples.
                                                # For example, for 3 columns, 'START', 'END', 'SIGNAL' - the first two are integers and the last a floating point decimal, you would have:
                                                # self.SQL = [('START','INT'),('END','INT'),('SIGNAL','REAL')]
                                                # Note that the type is directly given to the database, so you can go crazy if you want to use Postgres lists or something.

                                                # The data that will go into this table must be a list-of-lists, and like all stat modules, should be found under the same name as
                                                # the stat module itself. So say the module for the previous table layout was called "pileup", the data might look like:
                                                # pileup = [ [0,100,1], [101,200,3], .... [8382300,8382400,1] ]

        self.index = None                       # First of all, this parameter has no effect for linked stats, or json_stats. It only runs if you define the SQL table yourself.
                                                # You do not have to index your tables - sometimes it doesn't make sense to - but if you want to, here you would specify the SQL
                                                # to do it :) Essentially this is an SQL command run after the table has been created/populated.
                                                # If this is set to None, or self.index is missing entirely, no index will be made (but you can always do it manually later!)

        self.viz = ['catagorical',None]         # Another optional parameter for a stat. This one defines how to display data when visualized in the browser.
                                                # For linked stats, the default is to use self.SQL to determine if the data is catagorical ('TEXT') or continuous ('INT/REAL').
                                                # This can be over-ridden by specifying a value for self.viz here. 
                                                # The first element in the list is the name of the vizulization. The second element is the code to produce it.
                                                # If the code is not specified (second element has a value of None, False, or the list only has 1 element), but a name is
                                                # supplied, then the code to visualize the data will be pulled down from AC.GT. This is typically better when you want to share
                                                # your stat methods, and you dont want new versions of vizulization code to be effect your stats MD5 checksum. More importantly,
                                                # people dont typically trust external code in their browsers. Thus, you are encouraged to let us check your code and if its not
                                                # malicious, host it on AC.GT under a unique name. You can see all available visualizations on http://ac.gt/seqc/viz

                                                # However, if you want to impliment your own visualizations for yourself, or AC.GT hasn't gotten around to reviewing your code 
                                                # yet, then you can put in your own code here. This will be stored in the METHOD table of your database and can be queried by
                                                # other users if you make your database publicly availible (so you can share methods, data, and how to view that data!)
                                                # Users however will be warned that the vizulization code has not been checked by AC.GT, and they can accept (or not) the code
                                                # after reviewing it for themselves.

        self.dependancies = ['SEQ']             # Due to the way python interacts with hts-lib (the C code that actually gets data out of the BAM file), repeat requests
                                                # to something like "read.seq" (the DNA in the read's SEQ field) add up. If two stats thenfore both need SEQ to claculate
                                                # their statistic, we will be 2x slower processing the data than if we fetched SEQ once, then used it in both modules.
                                                # For this reason, we urge module developers to never interact with the underlying pysam/hts-python libraries directly if
                                                # they can help it. Instead, specify the stat name (in this case "SEQ" here. This way, the "SEQ" variable, will have it's 
                                                # value set and availible before this module is run. It is always better to use dependancies that roll your own, and due to the
                                                # way code is dynamically inlined, there's no penalty at all for using them. This means it can often be a good idea to break your
                                                # module up into a chain of modules all dependant on the previous modules in the chain. The TYPE module is a good example of this.
                                                # Dependancies is an optional parameter and can be excluded or set to None/False. If you do wish to provide one or more
                                                # Dependancies however, you must provide here a list of strings, where the string is the name of the module to run first.

        self.METHOD       = 'print "whoop!"'    # The code (as a string) that SeQC should execute for each read it processes in the SAM/BAM file is found in self.METHOD
                                                # The code can be many lines and do many things - the only criteria is that the final "result" for the stat (per read for 
                                                # linked stats, or at the very end for unlinked stats) must be stored under the name of the stat itself.
                                                # You can use indentations in any format you like. Such freedoms. What a time to be alive.

        if INFO['fileReader'] == 'sam':         # However, for some modules we need to know exactly how we are reading the file. This is particularly true for the built-in
          self.METHOD = 'QNAME = read[0]'       # functions where we cannot simply use other modules as dependancies (which is usually recommended).
        elif INFO['fileReader'] == 'pysam':     # We do this by passing the "fileReader" variable to the class when we initialize it, and then using that here to decide what
          self.METHOD = 'QNAME = read.qname'    # the method should be. Currently, values for fileReader can be 'sam', 'pysam' and 'htspython'. There may be more in the future.
        elif INFO['fileReader'] == 'htspython': # In fact the INFO my be used in the future to extend modules by setting extra paramters at runtime beyond fileReader.
          self.METHOD = 'QNAME = read.qname'    # If you decide to use multiple methods in your own stat, please make sure your stat always produces *exactly* same output!
        else: self.METHOD = None                # And if you use 'if' like we have here, make sure you end with an 'else: None'!

        self.before = None                      # Some modules require extra processing before reading the file. See the FLAG module as an example. This works exactly
                                                # like self.METHOD in that it should be a string. This will be run once and only once, no matter how many times your module is
                                                # called from the command line in multiple groups (because all groups are calculated simultaniously, with no redundancy)

        self.after  = None                      # Some modules require extra processing to be done after the file has been read. See the RNAME module as an example.
                                                # Unlike .before, the code here will run once per occurence of the stat on the command line. This is because .after is often 
                                                # used to edit/transform the data after it has all been collected. So where is the data?
                                                # For linked modules, data is stored under a variable called "table". table is the whole table of data before it gets put into 
                                                # the SQL database. It will contain as many columns as there were in the analysis group. To get to the data your stat made, 
                                                # use the "column" variable, which has been set to the id of your data's column in the row. eg. table[column][0]  would be the 
                                                # first value for the first row in the table that your module created.

                                                # IMPORTANT!! Modules that use self.after confuse people when they try to import them as a dependancy and they don't get back
                                                # something that looks like the final result they are used to - if you use .after, warn users about this in your module 
                                                # comments! Consider writing a version of your module that can be used as a dependancy if it needs to be, such as CHR is for 
                                                # RNAME.

addStat('STAT_NAME',[])                         # Finally, to initialize your stat with a name, and a list of MD5 sums of modules it is compatible with.
                                                # The stat name has to be short, and I recommend you don't use any symbols if you can help it. It MUST be the same as the result
                                                # that the stat makes per read/file.

                                                # The last thing, a list of MD5 sums of compatible modules, is something that is a bit confusing, but very important.
                                                # We want to keep our data consistent, such that we never end up comparing old stat results against new stat results, where
                                                # since the version change/update, the stat results are now somehow incompatible.
                                                # This is particularly important when two people from across the world create a stat called "depth", where the former doesn't 
                                                # include duplicate reads but the latter does. To get around this issue, all modules are MD5'd before going into the database's
                                                # METHODs table. When data is shared, so is the code to generate it. However, unlike a regular file MD5, only the data before 
                                                # the addStat() line is used in the hash. Everything else is deleted, and the stat name and list of MD5s are stored seperately.
                                                # This means changing the name of a stat, or setting dependancies, does not effect the MD5 of the stat itself. This allows two
                                                # module writers a way to mark each other's stat as compatible, or one of them can change their name with no fuss.

                                                # To get the MD5 sum of a module, use --md5 with it loaded as a .stat module, or look it up in the METHOD table of the 
                                                # resultant database.