# Load functions and packages
import helpers as hp
import pandas as pd
import os
from Bio import Entrez
from datetime import datetime, timedelta
import time
import numpy as np
import re
from urllib.error import HTTPError

issn_list = hp.process_csv()
hp.ask_email()

pubmed_id_df = hp.fetch_pubmed_ids_for_issns(issn_list)
pubmed_id_df.to_csv("data/pubmed_id_df.csv")
print("Saved PubMed IDs to data/pubmed_id_df.csv")

pubmed_id_df = pd.read_csv('data/pubmed_id_list.csv', usecols = ["PubMedID"], dtype = str)
# Post 5000 articles at a time then fetch batches of 10 iteratively
articles_batches = hp.fetch_articles(pubmed_id_df)
# Extract emails from each article
emails_df = hp.parse_author_emails(articles_batches)
# Save e-mails
emails_path = os.path.abspath("data/emails_df.csv")
emails_df.to_csv(emails_path)
print(f"Saved emails to {emails_path}")

hp.filter_emails(emails_df)
print("Filtered emails")
hp.create_master_df()
print("Created master dataframe")
