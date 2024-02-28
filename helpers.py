# Load packages
import pandas as pd
import os
from Bio import Entrez
from datetime import datetime, timedelta
import time
import numpy as np
import re
from urllib.error import HTTPError

# This function reads the CSV file and returns the top 10 journals based on the 5-year JIF for each category

def process_csv():
    current_dir = os.path.dirname(__file__)
    csv_path = os.path.join(current_dir, 'data', 'all_journals.csv')
    
    df = pd.read_csv(csv_path)
    
    df['category'] = df['category'].str.split(',')
    
    df = df.explode('category').reset_index(drop=True)

    grouped = df.groupby('category')
    top_jifs = grouped['JIF_5_year'].nlargest(10)
    
   # Reset the index of top_jifs to align with df's index
    top_jifs = top_jifs.reset_index(level=0, drop=True)
    
    # Retrieve the original rows for the top   10 JIFs in each category
    top_rows = df.loc[top_jifs.index]

    top_rows = top_rows.drop_duplicates(subset='ISSN').dropna(subset=['ISSN'])
    
    # Return the 'ISSN' values as a list
    return top_rows['ISSN'].tolist()

def ask_email():
    Entrez.email = input("Please enter your email address to use with requests: ")
    return None

def fetch_pubmed_ids_for_issns(issn_list):
    one_year_ago = (datetime.now() - timedelta(days=365)).strftime('%Y/%m/%d')
    pubmed_ids = []

    for issn in issn_list:
        start =  0
        max_results =  1000
        while True:
            try:
                # Search PubMed for articles with the given ISSN
                handle = Entrez.esearch(
                    db='pubmed',
                    term=f'{issn}[ta] AND {one_year_ago}[dp] :   2024/02/13[dp]',
                    retstart=start,
                    retmax=max_results
                )
                search_results = Entrez.read(handle)
                print(f'Current ISSN: {issn}')

                # Check if there are any results for the ISSN
                if not search_results['IdList']:
                    break

                # Add PubMed IDs to the list
                pubmed_ids.extend([(pubmed_id, issn) for pubmed_id in search_results['IdList']])

                
                # Update the start index for the next iteration
                start += max_results

            except Exception as e:
                print(f"An error occurred while searching for ISSN {issn}: {e}")
                break
        
        # Log progress
        print(f"Found {len(pubmed_ids)} articles so far")

    # Convert the list of tuples to a DataFrame
    pubmed_df = pd.DataFrame(pubmed_ids, columns=['PubMedID', 'ISSN'])

    # Remove rows with missing values
    pubmed_df = pubmed_df.dropna()

    # Remove rows with 0 values
    pubmed_df = pubmed_df[pubmed_df['PubMedID'] != 0]

    # Ensure PubMed IDs are strings
    pubmed_df['PubMedID'] = pubmed_df['PubMedID'].astype(str)

    # Print overall statistics

    print(f"Found {len(pubmed_df)} articles in total.")

    return pubmed_df

def fetch_articles(pubmed_df):
    error_encountered = False
    # Post 5000 ids at a time to avoid the 10k limit
    post_batch = 5000
    results = []
    batch_size = 10
    for i in range(0, len(pubmed_df), post_batch):
        search_results = Entrez.read(Entrez.epost("pubmed", id=",".join(pubmed_df['PubMedID'].astype(str)[i:i+post_batch])))
        print(f"Posting records {i} to {i+post_batch}...")
        for start in range(0, post_batch, batch_size):
            if error_encountered:
                print("You have likely fetched all posted articles. Stopping the fetching process.")
                break
            end = min(post_batch, start + batch_size)
            print(f"Fetching records {start} to {end}...")
            attempt = 1
            while attempt <= 3:
                try:
                    fetch_handle = Entrez.efetch(db="pubmed",
                                                 rettype="medline",
                                                 retmode="XML",
                                                 retstart=start,
                                                 retmax=batch_size,
                                                 webenv=search_results["WebEnv"],
                                                 query_key=search_results["QueryKey"])
                    data = Entrez.read(fetch_handle)
                    fetch_handle.close()
                    break
                except HTTPError as e:
                    # Check if the exception is an HTTPError with status code  400
                    if e.code == 400:
                        print(f"{e}")
                        error_encountered = True # Stop fetching if a  400 Bad Request error is encountered
                        attempt += 1
                        time.sleep(2)
                        continue
                except Exception as ex:
                    print(f"Attempt {attempt} failed with error: {ex}")
                    attempt +=1
                    time.sleep(2)
            if data:
                results.append(data)
    return results

def parse_author_emails(results):
    pattern = r'[a-z0-9\.\-+_]+@[a-z0-9\.\-+_]+\.[a-z]+'
    data_list = []  # Collect data in a list

    for data in results:
        # Get the PubmedArticle array
        pubmedarticles = data.get('PubmedArticle', [])
        # Get the PubmedBookArticle array
        pubmedbookarticles = data.get('PubmedBookArticle', [])

        # Iterate over each article
        for article in pubmedarticles:
            title = article['MedlineCitation']['Article']['ArticleTitle']
            pubmed_id = str(article['MedlineCitation']['PMID'])
            if 'AuthorList' not in article['MedlineCitation']['Article']:
                data_list.append({"pubmedid": pubmed_id, "last_name": np.NaN, "email": np.NaN, "message": 'No authors found.', "title": title})
                continue

            author_list = article['MedlineCitation']['Article']['AuthorList']
            for author in author_list:
                if 'LastName' not in author:
                    last_name = np.NaN
                else:
                    last_name = author['LastName']
                if 'AffiliationInfo' not in author:
                    data_list.append({"pubmedid": pubmed_id, "last_name": last_name, "email": np.NaN, "message": 'No affiliation found.', "title": title})
                    continue
                affiliation_list = author['AffiliationInfo']
                for affiliation in affiliation_list:
                    email_match = re.search(pattern, affiliation['Affiliation'])
                    if email_match:
                        data_list.append({"pubmedid": pubmed_id, "last_name": last_name, "email": email_match.group(), "message": 'None', "title": title})
                    else:
                                         data_list.append({"pubmedid": pubmed_id, "last_name": last_name, "email": np.NaN, "message": 'No email found', "title": title})

        # Iterate over each book article
        for bookarticle in pubmedbookarticles:
            pubmed_id = str(bookarticle['BookDocument']['PMID'])
            data_list.append({"pubmedid": pubmed_id, "last_name": np.NaN, "email": np.NaN, "message": 'PubmedBookArticle', "title": np.NaN})

        print(f"Processed {len(pubmedarticles)} PubmedArticle(s) and {len(pubmedbookarticles)} PubmedBookArticle(s).")

    df = pd.DataFrame(data_list)
    df['email'] = df['email'].replace('NaN', np.nan)
    unique_emails = df['email'].nunique()
    print(f"Found {unique_emails} unique email addresses.")
    # Ensure variable types
    df['pubmedid'] = df['pubmedid'].astype(str)
    df['last_name'] = df['last_name'].astype(str)
    df['email'] = df['email'].astype(str)
    df['message'] = df['message'].astype(str)
    df['title'] = df['title'].astype(str)
    return df

# Email filtering function

def filter_emails():
    emails_df = pd.read_csv('data/emails_df.csv')

    # Aggregate last names of authors
    aggregation = {
        'last_name': lambda x: ', '.join(x)
    }

    # Group by pubmedid, then aggregate
    df_grouped = emails_df.groupby(['pubmedid', 'email']).agg(aggregation).reset_index()

    # Filter rows with empty emails
    df_filtered = df_grouped[df_grouped['email'].notna()]

    # Drop the message column if you decide to keep it beforehand
    # df_filtered = df_filtered.drop(columns=['message'])

    # Save it
    df_filtered.to_csv("data/emails_filtered.csv")

def create_master_df():
    pubmed_id_df = pd.read_csv("data/pubmed_id_list.csv")
    pubmed_id_df = pubmed_id_df.rename(columns={'PubMedID' : "pubmedid", 'ISSN' : 'issn'})
    df_filtered = pd.read_csv("data/emails_filtered.csv")
    merged_df = pd.merge(df_filtered, pubmed_id_df, on='pubmedid', how='inner')
    journal_df = pd.read_csv('data/all_journals.csv')[["journal_name", "ISSN", "category"]]
    journal_df = journal_df.rename(columns={'ISSN' : 'issn'})
    merged_master = pd.merge(merged_df, journal_df, on='issn', how='left')
    merged_master = merged_master.loc[:, ~merged_master.columns.str.contains('Unnamed')]
    merged_master.to_csv("data/merged_master.csv", index=False)
