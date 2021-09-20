# BRENDA

## Donwloading database

Go to

```
https://www.brenda-enzymes.org/download_brenda.php#download
```

Register and download the file (this download the latest version).

Untar the file 

```
tar xvfz brenda_download.tar.gz
```


## Parsing

Go to the brenda directory
```
cd __YOUR_PATH__/brenda
```

Clone the external module
```
git clone https://github.com/alexandra-zaharia/BRENDA-Parser.git
```

Install the requirements (python >=3.6)

```
pip install -r requirements.txt
```

The parser contains a help
```
python parser_txt.py [-h] [--input INPUT_FILE] [--output OUTPUT_FILE]
```

The database can be parsed as
```
python parser_txt.py --input brenda_download.txt --output brenda_download.jsonl
```

The output file can be then uploaded to a CPS collection
