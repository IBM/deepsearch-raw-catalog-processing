# UNIPROT

## Donwloading database


Go to the uniprot directory

```
cd __YOUR_PATH__/uniprot
```

execute the script (this download the latest version)

```
./download_uniprot.x
```

## Parsing


Install the requirements (python >=3.6)

```
pip install -r requirements.txt
```

The parser contains a help
```
python parser.py [-h] [--input INPUT_FILE] [--output OUTPUT_FILE]
```

The database can be parsed as
```
python parser.py --input uniprot_sprot.xml.gz --output uniprot_sprot.jsonl
```

The output file can be then uploaded to a CPS collection
