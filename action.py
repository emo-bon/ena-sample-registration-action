import logging
import os
import pandas as pd
import re
import requests
import xml.etree.ElementTree as et
import yaml
from collections import namedtuple
from pathlib import Path
from uritemplate import expand
from urllib.error import HTTPError

GITHUB_WORKSPACE = Path(os.getenv("GITHUB_WORKSPACE", "/github/workspace"))

Endpoint = namedtuple("Endpoint", "username password dropbox query env")

ENDPOINT_PRO = Endpoint(
    os.getenv("USERNAME"),
    os.getenv("PASSWORD"),
    os.getenv("DROPBOX_PRO"),
    os.getenv("QUERY_PRO"),
    "PRO"
) 

ENDPOINT_DEV = Endpoint(
    os.getenv("USERNAME"),
    os.getenv("PASSWORD"),
    os.getenv("DROPBOX_DEV"),
    os.getenv("QUERY_DEV"),
    "DEV"
) 

def generate_ena_accession_numbers(smi, ref_code, endpoint):
    # TODO: replace github URLs with published rocrate urls in domain data.emobon.embrc.eu
    # TODO: when samples are triplized, get sample metadata from RDF graph rather than df/csv
    observatory = smi.split("_")[1]
    habitat = {"So": "sediment", "Wa": "water"}[smi.split("_")[2]]
    df_observatory = retrieve_metadata(
        "https://raw.githubusercontent.com/emo-bon/observatory-"
        f"{observatory.lower()}-crate/main/logsheets/{habitat}_observatory.csv",
        smi
    )
    df_sampling = retrieve_metadata(
        "https://raw.githubusercontent.com/emo-bon/observatory-"
        f"{observatory.lower()}-crate/main/logsheets/{habitat}_sampling.csv",
        smi,
        "source_mat_id"
    )
    df_governance = retrieve_metadata(
        "https://raw.githubusercontent.com/emo-bon/governance-data/main/observatories.csv",
        observatory,
        "EMOBON_observatory_id"
    )
    if not (df_observatory is None or df_sampling is None or df_governance is None):
        submission = get_submission_xml()
        sample = get_sample_xml(smi, ref_code, habitat, df_observatory, df_sampling)
        ean_sample = get_ean_from_ebi(smi, ref_code, sample, submission, endpoint)
        ean_project = df_governance["ENA_accession_number_project"].iloc[0]
        ean_umbrella = df_governance["ENA_accession_number_umbrella"].iloc[0]
        return ean_sample, ean_project, ean_umbrella

def retrieve_metadata(url, value, filter_column=None):
    try:
        df = pd.read_csv(url)
        if filter_column:
            df = df[df[filter_column] == value]
        assert df.shape[0] == 1
        return df
    except HTTPError as e:
        logger.error(f"HTTPError | Could not retrieve metadata for {value} ({url}) | {e}")
    except AssertionError as e:
        logger.error(f"AssertionError | Metadata size is not equal to 1 row for {value} ({url}) | {e}")


def get_submission_xml():
    submission = et.Element("SUBMISSION")
    actions = et.SubElement(submission, "ACTIONS")
    action = et.SubElement(actions, "ACTION")
    et.SubElement(action, "ADD")
    return submission


def get_sample_xml(smi, ref_code, habitat, df_observatory, df_sampling):
    sample_set = et.Element("SAMPLE_SET")
    sample = et.SubElement(sample_set, "SAMPLE")
    sample.attrib["alias"] = ref_code
    title = et.SubElement(sample, "TITLE")
    title.text = smi
    sample_name = et.SubElement(sample, "SAMPLE_NAME")
    taxon_id = et.SubElement(sample_name, "TAXON_ID")
    taxon_id.text = str(int(df_sampling["tax_id"].iloc[0])) # TODO: fix this in QC, should be int instead of float
    scientific_name = et.SubElement(sample_name, "SCIENTIFIC_NAME")
    scientific_name.text = df_sampling["scientific_name"].iloc[0]
    sample_attributes = et.SubElement(sample, "SAMPLE_ATTRIBUTES")
    if habitat == "sediment":  # ERC000021
        add_attribute(sample_attributes, "ENA-CHECKLIST", "ERC000021")
        add_attribute(sample_attributes, "elevation", str(0), "m")
    else:  # ERC000024
        add_attribute(sample_attributes, "ENA-CHECKLIST", "ERC000024")
    add_attribute(sample_attributes, "project name", str(df_observatory["project_name"].iloc[0]))
    add_attribute(sample_attributes, "collection date", str("2020-01-01")) # TODO: fix this in QC
    # add_attribute(sample_attributes, "collection date", str(df_sampling["collection_date"].iloc[0])) # TODO: fix this in QC
    add_attribute(sample_attributes, "geographic location (country and/or sea)", str(df_observatory["geo_loc_name"].iloc[0]))
    add_attribute(sample_attributes, "geographic location (latitude)", str(df_observatory["latitude"].iloc[0]), "DD")
    add_attribute(sample_attributes, "geographic location (longitude)", str(df_observatory["longitude"].iloc[0]), "DD")
    add_attribute(sample_attributes, "depth", str(df_sampling["depth"].iloc[0]), "m")
    add_attribute(sample_attributes, "broad-scale environmental context", str(df_observatory["env_broad_biome"].iloc[0]))
    add_attribute(sample_attributes, "local environmental context", str(df_observatory["env_local"].iloc[0]))
    add_attribute(sample_attributes, "environmental medium", str(df_sampling["env_material"].iloc[0]))
    return sample_set


def add_attribute(element, tag, value, units=None):
    sample_attribute = et.SubElement(element, "SAMPLE_ATTRIBUTE")
    t = et.SubElement(sample_attribute, "TAG")
    t.text = tag
    v = et.SubElement(sample_attribute, "VALUE")
    v.text = value
    if units:
        u = et.SubElement(sample_attribute, "UNITS")
        u.text = units


def get_ean_from_ebi(smi, ref_code, sample, submission, endpoint):
    files = {
        "SAMPLE": (
            "sample.xml",
            et.tostring(sample, xml_declaration=True, encoding="UTF-8")
        ),
        "SUBMISSION": (
            "submission.xml",
            et.tostring(submission, xml_declaration=True, encoding="UTF-8")
        ),
    }
    response = requests.post(
        endpoint.dropbox,
        auth=(endpoint.username, endpoint.password),
        files=files
    )
    assert response.headers["Content-Type"] == "application/xml"
    root = et.fromstring(response.content.decode())
    if root.attrib["success"] == "true":
        ean = root.find("SAMPLE").attrib["accession"] # TODO: was sample_ean
        # ext_id_ean = root.find("SAMPLE").find("EXT_ID").attrib["accession"]
        # submission_ean = root.find("SUBMISSION").attrib["accession"]
        return ean # TODO: decide which accession numbers to return
    else:
        try:
            if endpoint.env == "PRO":
                # TODO: test this block, as it is not available in DEV
                url = expand(endpoint.query, alias=f'"{ref_code}"')
                df = pd.read_csv(url, sep='\t')
                assert df.shape[0] == 1
                return df["sample_accession"].iloc[0]
            else:
                msg = root.find("MESSAGES").find("ERROR").text
                return re.search(r'"ERS\d{8}"', msg)[0][1:-1] # TODO: dangerous
        except Exception as e:
            logger.error(f"Exception | Could not register or retrieve ENA accession number for {smi} ({ref_code}) | {e}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    batches = [i for i in os.scandir(GITHUB_WORKSPACE / "shipment") if i.is_dir()]
    for batch in batches:
        meta = yaml.load(
            open(Path(batch.path) / f"metadata-{batch.name}.yml"),
            Loader=yaml.CLoader
        )
        if meta["ready_for_processing"] == True:
            if meta["production_deployment"] == True:
                endpoint = ENDPOINT_PRO
                ena_csv_path = Path(batch.path) / f"ena-accession-numbers-{batch.name}.csv"
            else:
                endpoint = ENDPOINT_DEV
                ena_csv_path = Path(batch.path) / f"ena-accession-numbers-{batch.name}.dev.csv"
            df = pd.read_csv(Path(batch.path) / f"list-of-samples-{batch.name}.csv")
            if ena_csv_path.exists():
                df_ena = pd.read_csv(ena_csv_path)
            else:
                df_ena = pd.DataFrame(
                    columns=[
                        "ref_code",
                        "source_material_id",
                        "ena_accession_number_sample",
                        "ena_accession_number_project",
                        "ena_accession_number_umbrella"
                    ]
                )
            for _, row in df.iterrows():
                ref_code = row["ref_code"]
                smi = row["source_material_id"]
                if not smi in df_ena["source_material_id"].to_list():
                    eans = generate_ena_accession_numbers(smi, ref_code, endpoint)
                    if eans and all(eans): # eans is not None and none of the elements in eans is None 
                        df_ena = df_ena.append(
                            {
                                "ref_code": ref_code,
                                "source_material_id": smi,
                                "ena_accession_number_sample": eans[0],
                                "ena_accession_number_project": eans[1],
                                "ena_accession_number_umbrella": eans[2],
                            },
                            ignore_index=True
                        )
                        df_ena.to_csv(ena_csv_path, index=False)
                        logger.info(f"ENA accession numbers generated for {smi} ({ref_code})")
                    else:
                        logger.info(f"No ENA accession numbers generated for {smi} ({ref_code})")
