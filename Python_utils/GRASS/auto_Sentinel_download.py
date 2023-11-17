'''
Auto-dowloading Sentinel data
Author: Songyan Zhu
Contact: szhu4@ed.ac.uk
---------Log of Changes-------------
Created: 2023-11-14
Updated: 2023-11-15
    |-> Make it operational
Updated: 2023-11-16
    |-> Only retrieve L2A product
    |-> Check file exisitence
    |-> Add error catch
    |-> Fix the '-1' remaining counting issue
'''

# Import credentials
import yaml
import requests
import argparse
import pandas as pd
import geopandas as gpd
from pathlib import Path
from creds import *

def load_yaml_config(p):
    with open(p, "r") as stream:
        try:
            yamlfile = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            assert(exc)
    return yamlfile

def get_access_token(username: str, password: str) -> str:
    data = {
        "client_id": "cdse-public",
        "username": username,
        "password": password,
        "grant_type": "password",
    }
    try:
        r = requests.post(
            "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token",
            data=data,
        )
        r.raise_for_status()
    except Exception as e:
        raise Exception(
            f"Access token creation failed. Reponse from the server was: {r.json()}"
        )
    return r.json()["access_token"]

def load_config(p):
    cfg = load_yaml_config(p)
    username = cfg['username']
    password = cfg['password']
    start_date = cfg['start_date']
    end_date = cfg['end_date']
    data_collection = cfg['data_collection']
    roi_str = cfg['roi_str']
    roi_file = cfg['roi_file']

    # ----------------------------------------------------------------------------------
    if roi_file:
        shp = gpd.read_file(roi_file)
        if shp.crs.to_epsg() != 4326: shp = shp.to_crs(4326)
        roi = shp.geometry.to_wkt()[0] + "'" # only takes the first row!
    elif roi_str:
        roi = roi_str
    else:
        raise Exception('Either roi_file or roi_str must exist in the configuration!')

    return username, password, start_date, end_date, data_collection, roi

def retrieve(p_config):
    username, password, start_date, end_date, data_collection, roi = load_config(p_config)

    json = requests.get(f"https://catalogue.dataspace.copernicus.eu/odata/v1/Products?$filter=Collection/Name eq '{data_collection}' and OData.CSC.Intersects(area=geography'SRID=4326;{roi}) and ContentDate/Start gt {start_date}T00:00:00.000Z and ContentDate/Start lt {end_date}T00:00:00.000Z").json()
    dfd = pd.DataFrame.from_dict(json['value'])
    dfd = dfd[dfd['S3Path'].str.contains('/L2A/')] # Only L2A products
    n_files = len(dfd)
    print(f'Beginning to retrieve {n_files} {data_collection} files...')
    for cnt in dfd.index:
        try:
            image_id = dfd.loc[cnt, 'Id']
            savename = data_collection + '-' + str(cnt).zfill(4) + '-' + dfd.loc[cnt, 'OriginDate'].split('T')[0].replace('-', '') + '.zip'
            if Path(savename).exists(): continue

            url = f"https://zipper.dataspace.copernicus.eu/odata/v1/Products({image_id})/$value"
            access_token = get_access_token(username, password)
            headers = {"Authorization": f"Bearer {access_token}"}

            session = requests.Session()
            session.headers.update(headers)
            response = session.get(url, headers=headers, stream=True)

            with open(savename, "wb") as file:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        file.write(chunk)
            print(f'{cnt + 1} done, {n_files - cnt + 1} remaining...')
            print('-' * 100)
        except Exception as e:
            print(e)
            with open('auto_Sentinel_download_logging.txt', 'a') as f:
                f.write(e + '|' + savename + '|' + image_id + 'n')
        
if __name__ == '__main__':
    # Example: python auto_Sentinel_download.py
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", nargs = "?", default = 'auto_Sentinel_download.yaml', type = str)
    args = parser.parse_args()
    p_config = args.p
    retrieve(p_config)