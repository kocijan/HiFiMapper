# read everything from the first file "sd_0001.fastq" into variable called "reads"
reads = []
with open("/home/mjkocijan/hifimapper/HiFiMapper/tests/benchmarks/data/sd_0001.fastq", "r") as f:
    for line in f:
        reads.append(line)

# read everything from the second file "sd_0001.info" into variable called "info"
info = []
with open("/home/mjkocijan/hifimapper/HiFiMapper/tests/benchmarks/data/sd_0001.info", "r") as f:
    for line in f:
        info.append(line)

j = 0
# go through reads line by line
for i in range(len(reads)):
    # if the line starts with "@S" and is shorter than 100 characters
    if reads[i].startswith("@S") and len(reads[i]) < 100:
        # instead of printing from the file "sd_0001.fastq", print from the file "sd_0001.info" line j
        # the line j might like like "27704395 10593 +" and we want to print
        # "@CP068265.2_0_TARGET_SIZE_10593_STRAIN_+_POSITION_27704395"
        try:
            print("@CP068265.2_0_TARGET_SIZE_" + info[j].split()[1] + "_STRAIN_" + info[j].split()[2] + "_POSITION_" + info[j].split()[0])
            pass
        except IndexError:
            exit(0)
        j += 1
    else:
        print(reads[i], end="")
        pass
