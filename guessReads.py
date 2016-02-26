'''
This little function used to go in place of MD5andCount().
It guesses the number of reads in a BAM/SAM file based on how many reads in the first X bytes.
It was used to produce accurate-ish status bars when processing the BAM/SAM files, however now SeQC
gets an accurate read count 'for free' when MD5ing the file, so its not used.
Left here for anyone else looking to get read-count estimates. Which is probably no one.
'''

# Estimate reads for SAM and BAM files
def guessReads(fileName,samtools):
  DEVNULL = open(os.devnull, 'wb')
  bytesUsed = 10000000
  #sizeInBytes = int(subprocess.check_output("ls -ln '" + str(fileName) + "' | awk '{print $5}'", stderr=DEVNULL, shell=True))
  sizeInBytes = os.path.getsize(fileName)
  if bamCheck(fileName):
      if bytesUsed*2 < sizeInBytes: # because if we can read the whole file in twice the time as it took to sample it, we might as well just read the whole file!
          readsIn100MB = int(subprocess.check_output('head -c ' + str(bytesUsed) + ' "' + str(fileName) + '" | "' + samtools + '" view - | wc -l', stderr=DEVNULL, shell=True))
          totalReads = (readsIn100MB/float(bytesUsed))*sizeInBytes * 1.1 # (1.1 adds 10% to read estimate so we are always conservative)
      else:
          totalReads = int(subprocess.check_output('"' + samtools + '" view -c "' + str(fileName) + '"', stderr=DEVNULL, shell=True))
  else:
      bytesIn10000Lines = int(subprocess.check_output('head -11000 "' + str(fileName) + '" | tail -10000 | wc -c', stderr=DEVNULL, shell=True))
      if bytesIn10000Lines == 0: return 0
      totalReads = ( sizeInBytes/float(bytesIn10000Lines) )*10000
      totalReads = int(totalReads * 1.05) # becase we count bytes in 10000 lines (and no header), estimate is more accurate than BAM
  return int(totalReads)
totalReads = guessReads(inputFile,args.samtools)

'''
Definitely no one.
'''
