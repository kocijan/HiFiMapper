import os
import subprocess
import shutil


WORK_DIRECTORY = os.getcwd()
SIMULATOR_DIRECTORY = f"{WORK_DIRECTORY}/simulator"
MAPPER_DIRECTORY = f"{WORK_DIRECTORY}/src"
EXE_DIRECTORY = f"{WORK_DIRECTORY}/exe"


# Make exe directory
os.makedirs(EXE_DIRECTORY, exist_ok=True)

os.system(f"g++ -std=c++17 -fopenmp -I ./vendor/bioparser/include -g -pthread {MAPPER_DIRECTORY}/main.cpp -lm -lz -o HiFiMapper")
os.system(f"g++ -std=c++17 -fopenmp -I ./vendor/biosoup/include -I ./vendor/bioparser/include -g -pthread {SIMULATOR_DIRECTORY}/reads_generator.cpp -lm -lz -o Readsgen")
os.system(f"g++ -std=c++17 -fopenmp -g {SIMULATOR_DIRECTORY}/reference_generator.cpp -lm -lz -o Refgen")

shutil.copy2("HiFiMapper", EXE_DIRECTORY)
shutil.copy2("Refgen", EXE_DIRECTORY)
shutil.copy2("Readsgen", EXE_DIRECTORY)

os.remove("HiFiMapper")
os.remove("Refgen")
os.remove("Readsgen")
