# diplomski

```
python3 build.py
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

Run experiments in tests/benchmarks:
```
python3.8 benchmark[1...].py
```

