#!/bin/bash

FTP=ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/
TARGET_DIR=./


wget -P $TARGET_DIR $FTP/uniprot_sprot.xml.gz


