# diplomski

```
git clone https://github.com/SuzanaPratljacic/HiFiMapper.git && cd HiFiMapper && mkdir build && cd build
cmake -DBUILD_SIMULATOR=True .. && make

```

To run experiments, first install minimap and Winnowmap:

```
git clone https://github.com/lh3/minimap2
cd minimap2 && make
```

```
git clone https://github.com/marbl/Winnowmap.git
cd Winnowmap
make -j8
```

Install python requirements:
```
pip install requirements.txt
```

Preperae the data

cd tests/benchmarks
mkdir data && cd data
Download [data](https://drive.google.com/drive/folders/1agwq4nsYo649RzUhUuCPIfTOMc4C0dfO?usp=sharing).

Run experiments in tests/benchmarks:
```
cd tests/benchmarks 
python3.8 benchmark[1...].py
```

