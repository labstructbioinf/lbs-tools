'''generate embeddings from sequence'''
import re
import os
import argparse
import pandas as pd
from tqdm import tqdm
from transformers import T5Tokenizer, T5EncoderModel
import torch

regex_aa = re.compile(r"[UZOB]")
EMBEDDER = 'Rostlab/prot_t5_xl_half_uniref50-enc'
BATCH_SIZE = 64

parser = argparse.ArgumentParser(description =  
    """
    Embedding script create embeddings from sequences via prot_t5_xl_half_uniref50
    by default `seq` column in used as embedder input. Records are stored
    as list maintaining dataframe order. 
    In python load via: 
    >>> import torch
    >>> torch.load(..)
    or 
    >>> import pickle
    >>> with open(.., 'rb') as f:
    >>>    embs = pickle.load(f)
    
    example use:
        python embeddings.py df.csv df.pt -cname seqfull
    """,
    formatter_class=argparse.RawDescriptionHelpFormatter
    )
parser.add_argument('input', help='csv/pickle (.csv or .p) with `seq` column',
                    type=str)
parser.add_argument('output', help='resulting list of embeddings',
                    type=str)
parser.add_argument('-cname', help='custom sequence column name',
                     dest='cname', type=str, default='')
parser.add_argument('-r', '-head', help='number of rows from begining to use',
                    dest='head', type=int, default=0)
args = parser.parse_args()
if args.input.endswith('csv'):
    df = pd.read_csv(args.input)
elif args.input.endswith('.p'):
    df = pd.read_pickle(args.input)
else:
    raise FileNotFoundError(f'invalid input infile extension {args.input}')

out_basedir = os.path.dirname(args.output)
if out_basedir == '':
    pass
else:
    if not os.path.isdir(out_basedir):
        raise FileNotFoundError(f'output directory is invalid: {out_basedir}')


if args.cname != '':
    if args.cname not in df.columns:
        raise KeyError(f'no column: {args.cname} available in file: {args.input}')
    else:
        print(f'using column: {args.cname}')
        if 'seq' in df.columns:
            df.drop(columns=['seq'], inplace=True)
        df.rename(columns={args.cname: 'seq'}, inplace=True)

if args.head > 0:
    df = df.head(args.head)
df.reset_index(inplace=True)

print('loading models')
tokenizer = T5Tokenizer.from_pretrained(EMBEDDER, do_lower_case=False)
model = T5EncoderModel.from_pretrained(EMBEDDER, torch_dtype=torch.float32)
num_records = df.shape[0]
residues = num_records % BATCH_SIZE
num_batches = int(num_records/BATCH_SIZE)
if residues > 0:
    num_batches += 1
print('num seq:', num_records)
print('num batches:', num_batches)

embedding_stack = list()
for batch_id in tqdm(range(num_batches)):
    seqlist = []
    lenlist = []
    for i, (idx, row) in enumerate(df.iterrows()):
        # use only current batch sequences
        if batch_id*BATCH_SIZE <= i < (batch_id + 1)*BATCH_SIZE:
            sequence = row.seq
            sequence = regex_aa.sub("X", row.seq)
            lenlist.append(len(sequence))
            sequence = " ".join(list(sequence))
            seqlist.append(sequence)
        
    ids = tokenizer.batch_encode_plus(seqlist, add_special_tokens=True, padding="longest")
    input_ids = torch.tensor(ids['input_ids'])
    attention_mask = torch.tensor(ids['attention_mask'])
    
    with torch.no_grad():
        embeddings = model(input_ids=input_ids, attention_mask=attention_mask)
        embeddings = embeddings[0].float().cpu()
    # remove sequence padding
    num_batch_embeddings = len(embeddings)
    assert num_batch_embeddings == len(seqlist)
    for i in range(num_batch_embeddings):
        seq_len = lenlist[i]
        emb = embeddings[i]
        if emb.shape[1] < seq_len:
            raise KeyError('sequence is longer then embedding')
        emb_no_padding = emb[:seq_len]
        embedding_stack.append(emb_no_padding)

torch.save(embedding_stack, args.output)