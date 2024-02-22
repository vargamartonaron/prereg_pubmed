# Load packages
import pandas as pd
import os
from Bio import Entrez
from datetime import datetime, timedelta
import time
import numpy as np
import re

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

def fetch_pubmed_ids_for_issns(issn_list, email='martonaronvarga@gmail.com'):
    Entrez.email = email
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

def fetch_articles(pubmed_df, email='martonaronvarga@gmail.com'):
    Entrez.email = email
    # Post 5000 ids at a time to avoid the 10k limit
    post_batch = 5000
    results = []
    batch_size = 10
    for i in range(0, len(pubmed_df), post_batch):
        search_results = Entrez.read(Entrez.epost("pubmed", id=",".join(pubmed_df[i:i+post_batch])))
        print(f"Posting records {i} to {i+post_batch}...")
        for start in range(0, post_batch, batch_size):
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
                except Exception as e:
                    print(f"Attempt {attempt} failed with error: {e}")
                    attempt += 1
                    time.sleep(2)
            if data:
                results.append(data)
    return results

def parse_author_emails(results):
    pattern = r'[a-z0-9\.\-+_]+@[a-z0-9\.\-+_]+\.[a-z]+'
    df = pd.DataFrame(columns=["pubmedid", "email", "message"])

    for data in results:  
        # Get the PubmedArticle array
        pubmedarticles = data.get('PubmedArticle', [])
        # Get the PubmedBookArticle array
        pubmedbookarticles = data.get('PubmedBookArticle', [])

        # Iterate over each article
        for article in pubmedarticles:
            pubmed_id = str(article['MedlineCitation']['PMID'])
            if 'AuthorList' not in article['MedlineCitation']['Article']:
                df = pd.concat([df, pd.DataFrame([{"pubmedid": pubmed_id, "last_name": np.NaN, "email": np.NaN, "message": 'No authors found.'}])])
                continue

            author_list = article['MedlineCitation']['Article']['AuthorList']
            for author in author_list:
                if 'LastName' not in author:
                    last_name = np.NaN
                else:
                    last_name = author['LastName']
                if 'AffiliationInfo' not in author:
                    df = pd.concat([df, pd.DataFrame([{"pubmedid": pubmed_id, "last_name": last_name, "email": np.NaN, "message": 'No affiliation found.'}])])
                    continue
                affiliation_list = author['AffiliationInfo']
                for affiliation in affiliation_list:
                    email_match = re.search(pattern, affiliation['Affiliation'])
                    if email_match:
                        df = pd.concat([df, pd.DataFrame([{"pubmedid": pubmed_id, "last_name": last_name, "email": email_match.group(), "message": 'None'}])])
                    else:
                        df = pd.concat([df, pd.DataFrame([{"pubmedid": pubmed_id, "last_name": last_name, "email": np.NaN, "message": 'No email found'}])])

        # Iterate over each book article
        for bookarticle in pubmedbookarticles:
            pubmed_id = str(bookarticle['BookDocument']['PMID'])
            df = pd.concat([df, pd.DataFrame([{"pubmedid": pubmed_id, "last_name": np.NaN, "email": np.NaN, "message": 'PubmedBookArticle'}])])

    df['email'] = df['email'].replace('NaN', np.nan)
    unique_emails = df['email'].nunique()
    print(f"Found {unique_emails} unique email addresses.")
    return df

