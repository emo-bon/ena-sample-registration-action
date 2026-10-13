cimport logging
import os
import pandas as pd
import requests
import time
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
    "https://www.ebi.ac.uk/ena/submit/drop-box/submit/",
    "https://www.ebi.ac.uk/ena/portal/api/search?result=sample&query=(sample_alias={alias})",
    "PRO"
) 

ENDPOINT_DEV = Endpoint(
    os.getenv("USERNAME"),
    os.getenv("PASSWORD"),
    "https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/",
    "https://wwwdev.ebi.ac.uk/ena/portal/api/search?result=sample&query=(sample_alias={alias})",
    "DEV"
) 

ERC000021 = pd.read_csv("https://raw.githubusercontent.com/emo-bon/sequencing-profile/main/checklist-translations/ERC000021.csv")
ERC000021 = ERC000021[ERC000021["include"] == "Y"]

ERC000024 = pd.read_csv("https://raw.githubusercontent.com/emo-bon/sequencing-profile/main/checklist-translations/ERC000024.csv")
ERC000024 = ERC000024[ERC000024["include"] == "Y"]


def generate_ena_accession_numbers(smi, ref_code, df_run_info, endpoint):
    # TODO: replace github URLs with published rocrate urls in domain data.emobon.embrc.eu
    # TODO: when samples are triplized, get sample metadata from RDF graph rather than df/csv
    observatory = smi.split("_")[1]
    habitat = {"So": "sediment", "Wa": "water"}[smi.split("_")[2]]
    df_sequencing = df_run_info[df_run_info["source_mat_id"] == smi]
    assert len(df_sequencing) == 1

    df_measured = retrieve_metadata(
        "https://raw.githubusercontent.com/emo-bon/observatory-"
        f"{observatory.lower()}-crate/refs/heads/main/logsheets/transformed/{habitat}_measured.csv",
        smi,
        "source_mat_id"
    )

    df_observatory = retrieve_metadata(
        "https://raw.githubusercontent.com/emo-bon/observatory-"
        f"{observatory.lower()}-crate/refs/heads/main/logsheets/transformed/{habitat}_observatory.csv",
        smi
    )

    df_sampling = retrieve_metadata(
        "https://raw.githubusercontent.com/emo-bon/observatory-"
        f"{observatory.lower()}-crate/refs/heads/main/logsheets/transformed/{habitat}_sampling.csv",
        smi,
        "source_mat_id"
    )

    df_governance = retrieve_metadata(
        "https://raw.githubusercontent.com/emo-bon/governance-data/main/observatories.csv",
        observatory,
        "EMOBON_observatory_id"
    )

    if not (df_measured is None or df_observatory is None or df_sampling is None or df_governance is None):
        submission = get_submission_xml()
        sample = get_sample_xml(smi, ref_code, habitat, df_measured, df_observatory, df_sampling, df_sequencing)
        ean_sample, an_biosamples = get_ean_from_ebi(smi, ref_code, sample, submission, endpoint)
        ean_project = df_governance["ENA_accession_number_project"].iloc[0]
        ean_umbrella = df_governance["ENA_accession_number_umbrella"].iloc[0]
        return ean_sample, ean_project, ean_umbrella, an_biosamples


def retrieve_metadata(url, value, filter_column=None):
    try:
        df = pd.read_csv(url, dtype=object, keep_default_na=False)
        if filter_column:
            df = df[df[filter_column] == value]
        assert len(df) == 1
        return df
    except HTTPError as e:
        logger.error(f"HTTPError | Could not retrieve metadata for {value} ({url}) | {e}")
    except AssertionError as e:
        logger.error(f"AssertionError | Metadata size is not equal to 1 row for {value} ({url}) | {e}")


def get_submission_xml():
    submission = et.Element("SUBMISSION")
    actions = et.SubElement(submission, "ACTIONS")
    action = et.SubElement(actions, "ACTION")
    et.SubElement(action, "ADD")  # TODO UPDATE?
    return submission


def get_sample_xml(smi, ref_code, habitat, df_measured, df_observatory, df_sampling, df_sequencing):
    df_sampling = df_sampling.copy()
    df_sampling["elevation"] = "0"
    sample_set = et.Element("SAMPLE_SET")
    sample = et.SubElement(sample_set, "SAMPLE")
    sample.attrib["alias"] = ref_code
    title = et.SubElement(sample, "TITLE")
    title.text = smi
    sample_name = et.SubElement(sample, "SAMPLE_NAME")
    taxon_id = et.SubElement(sample_name, "TAXON_ID")
    taxon_id.text = str(df_sampling["tax_id"].iloc[0].split("=")[1])
    scientific_name = et.SubElement(sample_name, "SCIENTIFIC_NAME")
    scientific_name.text = df_sampling["scientific_name"].iloc[0]
    sample_description = et.SubElement(sample, "DESCRIPTION")
    sample_description.text = df_sampling["samp_description"].iloc[0]
    sample_attributes = et.SubElement(sample, "SAMPLE_ATTRIBUTES")
    
    if habitat == "sediment":  # ERC000021
        add_attribute(sample_attributes, "ENA-CHECKLIST", "ERC000021")
        df_checklist = ERC000021
    else:  # ERC000024
        add_attribute(sample_attributes, "ENA-CHECKLIST", "ERC000024")
        df_checklist = ERC000024

    sheet2df = {
        "measured": df_measured,
        "observatory": df_observatory,
        "sampling": df_sampling,
        "sequencing": df_sequencing
    }

    for _, row in df_checklist.iterrows():
        tag = row["ENA_term"]
        value = str(sheet2df[row["EMO_BON_sheet"]][row["EMO_BON_term"]].iloc[0])
        units = row["units"]
        if pd.isna(units): units = None
        if tag == "sequence quality check": value = "manual"  # TODO get this sorted out in QC
        if tag == "tidal stage" and value == "ebb_tide": value = "low"  # TODO get this sorted out in QC
        if tag == "tidal stage" and value == "flood_tide": value = "high"  # TODO get this sorted out in QC
        if tag == "tidal stage" and value == "low_tide": value = "low"  # TODO get this sorted out in QC
        if tag == "tidal stage" and value == "high_tide": value = "high"  # TODO get this sorted out in QC
        if tag == "tidal stage" and value == "no_tide": value = ""  # TODO get this sorted out in QC
        if not value in ("", "NA", "nan"):
            add_attribute(sample_attributes, tag, value, units)

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
    time.sleep(0.1)
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
        # Save XML files locally for inspection
    xml_dir = GITHUB_WORKSPACE / "xml-output"
    xml_dir.mkdir(parents=True, exist_ok=True)

    sample_xml_path = xml_dir / f"{ref_code}_sample.xml"
    submission_xml_path = xml_dir / f"{ref_code}_submission.xml"

    et.ElementTree(sample).write(
        sample_xml_path,
        encoding="UTF-8",
        xml_declaration=True
    )

    et.ElementTree(submission).write(
        submission_xml_path,
        encoding="UTF-8",
        xml_declaration=True
    )

    logger.info(f"Saved XML files: {sample_xml_path}, {submission_xml_path}")

    try:
        response = requests.post(
            endpoint.dropbox,
            auth=(endpoint.username, endpoint.password),
            files=files,
            timeout=60
        )

        logger.info(f"ENA HTTP Status: {response.status_code}")
        logger.info(f"ENA Response Headers: {dict(response.headers)}")
        logger.info(f"ENA Response Body: {response.text[:3000]}")

        # Save the ENA receipt for inspection
        receipt_path = xml_dir / f"{ref_code}_receipt.xml"
        receipt_path.write_bytes(response.content)

        response.raise_for_status()

        # Validate XML response
        root = et.fromstring(response.content)

    except requests.exceptions.RequestException as e:
        logger.error(f"ENA HTTP error for {smi}: {e}")
        return None, None

    except et.ParseError as e:
        logger.error(f"ENA returned invalid XML for {smi}: {e}")
        return None, None

    if root.attrib["success"] == "true":
        sample_ean = root.find("SAMPLE").attrib["accession"]
        ext_id_ean = root.find("SAMPLE").find("EXT_ID").attrib["accession"]  # TODO if type="biosample"
        return sample_ean, ext_id_ean
    else:
        try:
            url = expand(endpoint.query, alias=f'"{ref_code}"')
            df = pd.read_csv(url, sep='\t')
            assert len(df) == 1
            return df["sample_accession"].iloc[0], "not found"
        except Exception as e:
            logger.error(f"Exception | Could not register or retrieve ENA accession number for {smi} ({ref_code}) | {e}")
            return None, None


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    batches = [i for i in os.scandir(GITHUB_WORKSPACE / "shipment") if i.is_dir()]
    for batch in batches:
        properties_path = Path(batch.path) / f"properties-{batch.name}.yml"
        if properties_path.exists():
            properties = yaml.load(
                open(properties_path),
                Loader=yaml.CLoader
            )
            if properties["ready_for_processing"] == True:
                # Select target genes
                if properties["target_gene"] == "metag":
                    TARGET_GENES = ["metag"]

                elif properties["target_gene"] == "metab":
                    TARGET_GENES = ["COI", "18S rRNA"]

                else:
                    logger.warning(
                        f"Unknown target_gene: {properties['target_gene']}"
                    )
                    continue

                if properties["production_deployment"] == True:
                    endpoint = ENDPOINT_PRO
                    ena_csv_path = Path(batch.path) / f"ena-accession-numbers-{batch.name}.csv"
		    
                else:
                    endpoint = ENDPOINT_DEV
                    ena_csv_path = Path(batch.path) / f"ena-accession-numbers-{batch.name}.dev.csv"
                df_run_info = pd.read_csv(Path(batch.path) / f"run-information-{batch.name}.csv")

                df_run_info = df_run_info[
                    df_run_info["target_gene"]
                    .astype(str)
                    .str.strip()
                    .str.lower()
                    .isin([gene.lower() for gene in TARGET_GENES])
                ].copy()

                logger.info(
                    f"Selected {len(df_run_info)} sequencing records "
                    f"for target genes: {TARGET_GENES}"
                )

                if ena_csv_path.exists():
                    df_ena = pd.read_csv(ena_csv_path)
                else:
                    df_ena = pd.DataFrame(
                        columns=[
                            "source_mat_id",
                            "ref_code",
                            "ref_code_seq",
                            "ena_accession_number_sample",
                            "ena_accession_number_project",
                            "ena_accession_number_umbrella",
                            "biosamples_accession_number"
                        ]
                    )
                for _, row in df_run_info.iterrows():
                    smi = row["source_mat_id"]
                    ref_code = row["ref_code"]
                    ref_code_seq = row["ref_code_seq"]
                    if (not smi in df_ena["source_mat_id"].to_list()) and ("BPNS" in smi): # TODO: remove BPNS filter
                        eans = generate_ena_accession_numbers(smi, ref_code, df_run_info, endpoint)
                        if eans and all(eans): # eans is not None and none of the elements in eans is None 
                            df_ena = df_ena.append(
                                {
                                    "source_mat_id": smi,
                                    "ref_code": ref_code,
                                    "ref_code_seq": ref_code_seq,
                                    "ena_accession_number_sample": eans[0],
                                    "ena_accession_number_project": eans[1],
                                    "ena_accession_number_umbrella": eans[2],
                                    "biosamples_accession_number": eans[3]
                                },
                                ignore_index=True
                            )
                            df_ena.to_csv(ena_csv_path, index=False)
                            logger.info(f"ENA accession numbers generated for {smi} ({ref_code})")
                        else:
                            logger.info(f"No ENA accession numbers generated for {smi} ({ref_code})")
