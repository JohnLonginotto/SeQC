# SeQC

_SeQC is still very much in the pre-Alpha stage of development. A lot has changed over the last year since the videos you may have seen on vimeo (http://vimeo.com/123508180), and i'm currently in the process of updating the website code to match the new database structure. Fortunately the database-creation side of things is fairly locked-in._

_There are a lot of new ideas for the website, including sharing QC stats/visualizations with others over the internet, bringing online all the ENCODE/Roadmap/etc QC data for comparison with your own, and a online modules directory to find trusted modules which do whatever you need them to do. That is all currently in development (or rather, will start up again after the PhD in about 2 months). OK - back to the blurb:_

SeQC is a "quality-control analysis" tool.

Most QC tools return reports, where the input data has either passed or failed the tool's QC tests.
Unfortunately, with high-throughput sequencing being the hammer for all nails that it is, there really is no way to write generic tests for quality, as programs such as FASTQC demonstrate. Many people will apply FASTQC to their data, and if it passes all the tests, assume that everything is OK. Alternatively, FASTQC or similar tools might fail your data on a specific test, but without anything to compare that failure too, the test itself becomes irrelevent and ignored. 

> "So what if my GC% distribution is out a bit - maybe thats biological?"

So SeQC is not like FASTQC or similar such reporting tools. Rather, it allows you to investigate the nature of sequencing data at scale, detect subtle inconsistencies by comparing your data with others, and ultimately really get to grips with the composition of your sequencing. It is also incredibly easy to use, and if you know some Python, incredibly easy to extend upon.

It also allows you to centralize your QC reports into a single database - either SQLite or Postegres - and even share/compare that data with others over the internet. How does your data compare to the ENCODE data sets? SeQC will tell you. Why does one of your replicates deviate so much from the other samples? SeQC will give you the tools to find out.

UPDATE 13th November 2016:
SeQC will be getting a big update soon, notably:

 1. pybam will be used instead of pysam/htspython, because its faster, 100% pypy complient, and works in the same way SeQC works (reading through the whole file once). This will also massively simplify the code, particularly for the BAM header, and SeQC will have 0 dependencies.
 2. Modules will go from being python classes to json objects.
 3. Modules will be able to have parameters (either required or optional), which will directly effect how argparse works so these parameters become 'native' when the module is loaded. For example ./SeQC --analysis GTF --GTF_FILE ./path/to/file.gtf
 4. The output SQL databases made by SeQC are becoming their own project - SeQL databases. SeQC will drop all the javascript code as a result, but still contain the code to make/add data to a SeQL database. All the other projects (Signl, ACGTrie, BAM+, etc) will also switch over to SeQL output, so literally everything will be in the same self-describing SQL format, with vizulization modules built-in.
 5. The master process will use select/fcntl as per the log.bio project to read subprocesses output rather than the ghetto ping.pong() method currently used. This is good because it means if a subprocess hangs or has a lengthy .before routine, it doesn't pause the status updating of all subprocesses. I anticipate .before and .after being more heavily used in the future, so this is important.
 6. SeQC will go up on ac.gt finally, i'll spend a full week on just documentation, and then try and publish it in the journal of bioinformatics with all the project contributors getting authorship.
