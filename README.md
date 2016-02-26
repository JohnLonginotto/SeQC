# SeQC

SeQC is a "quality-control analysis" tool.

Most QC tools return reports, where the input data has either passed or failed the tool's QC tests.
Unfortunately, with high-throughput sequencing being the hammer for all nails that it is, there really is no way to write generic tests for quality, as programs such as FASTQC demonstrate. Many people will apply FASTQC to their data, and if it passes all the tests, assume that everything is OK. Alternatively, FASTQC or similar tools might fail your data on a specific test, but without anything to compare that failure too, the test itself becomes irrelevent and ignored. 

> "So what if my GC% distribution is out a bit - maybe thats biological?"

So SeQC is not like FASTQC or similar such reporting tools. Rather, it allows you to investigate the nature of sequencing data at scale, detect subtle inconsistencies by comparing your data with others, and ultimately really get to grips with the composition of your sequencing. It is also incredibly easy to use, and if you know some Python, incredibly easy to extend upon.

It also allows you to centralize your QC reports into a single database - either SQLite or Postegres - and even share/compare that data with others over the internet. How does your data compare to the ENCODE data sets? SeQC will tell you. Why does one of your replicates deviate so much from the other samples? SeQC will give you the tools to find out.
