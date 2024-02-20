import helpers as hp
import os
import pandas as pd

out_path = os.path.abspath("data/pubmed_id_list.csv")

pubmed_id_df = pd.read_csv(out_path)
pubmed_id_df = pubmed_id_df['PubMedID'].map(str)
articles_batches = hp.fetch_articles(pubmed_id_df)
emails_df = hp.parse_author_emails(articles_batches)
emails_path = os.path.abspath("data/emails_df.csv")
emails_df.to_csv(emails_path)
