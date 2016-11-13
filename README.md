# SeQC

_SeQC is still very much in the Alpha stage of development. A lot has changed over the last year since the videos you may have seen on vimeo (http://vimeo.com/123508180), so it is advised to not get too familar with the SeQC code until after my PhD has finished and I can actually work on this and other projects full-time._

SeQC is a "quality-control analysis" tool, but it is not like FASTQC or similar. Written entirely in Python, SeQC is designed to help you to investigate the nature of sequencing data at scale, get to grips with your raw data, and write your own modules to do custom analysis with very little extra code needed.

It differs from other tools by outputing data to either an SQLite or Postegres SQL database (SeQL format). Data in SeQL format can be vizulised through the SeQL webservice, which comes with a number of generic plugins to look at quantitative data. All of your QC data from your project, your lab, or even your whole institute, can be loaded into a single SeQL database, so all your QC data is in one place.

The SeQL webservice will also allow you to communicate with other SeQL databases (to compare data) and offers you the ability to host your data for others to see - either publicly or privately. For SeQC this means you can compare your mapping rates to other publications, look for odd trends in your data, etc etc.

Once the code has been finalized and we have settled on the first Beta release, a full guide to using and extending SeQC will be published :) 

UPDATE 13th November 2016:
SeQC will be getting a big update soon, notably:

 1. pybam will be used instead of pysam/htspython, because its faster, 100% pypy complient, and works in the same way SeQC works (reading through the whole file once). This will also massively simplify the code, particularly for the BAM header, and SeQC will have 0 dependencies.
 2. Modules will go from being python classes to json objects.
 3. Modules will be able to have parameters (either required or optional), which will directly effect how argparse works so these parameters become 'native' when the module is loaded. For example ./SeQC --analysis GTF --GTF_FILE ./path/to/file.gtf
 4. The output SQL databases made by SeQC are becoming their own project - SeQL databases. SeQC will drop all the javascript code as a result, but still contain the code to make/add data to a SeQL database. All the other projects (Signl, ACGTrie, BAM+, etc) will also switch over to SeQL output, so literally everything will be in the same self-describing SQL format, with vizulization modules built-in.
 5. The master process will use select/fcntl as per the log.bio project to read subprocesses output rather than the ghetto ping.pong() method currently used. This is good because it means if a subprocess hangs or has a lengthy .before routine, it doesn't pause the status updating of all subprocesses. I anticipate .before and .after being more heavily used in the future, so this is important.
 6. SeQC will go up on ac.gt finally, i'll spend a full week on just documentation, and then try and publish it in the journal of bioinformatics with all the project contributors getting authorship.
