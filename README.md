The action that will register sample accessions in ENA. Input to this comes from 
* The files in the [batches](https://github.com/emo-bon/sequencing-data/tree/main/shipment) folders
* The [checklist translations](https://github.com/emo-bon/sequencing-profile)
* Project and umbrella ENA accession numbers from [https://github.com/emo-bon/governance-data/blob/main/observatories.csv](https://github.com/emo-bon/governance-data/blob/main/observatories.csv)
* Transformed logsheets found in each separate observatory crate

The action creates XML files which are then uploaded to ENA, returning sample accession numbers, which are then put back into the relevant files in the [batches](https://github.com/emo-bon/sequencing-data/tree/main/shipment) folders.
The code that sits behind this is called action.py

`LOCAL RUN`
Build docker: 

````sh
sudo docker build -t emobon-ena-register .
````
Prepare ENV:
````sh
export GITHUB_WORKSPACE="$(pwd)/ena-test"
export USERNAME="XXXXXX"
export PASSWORD="XXXXXX"
````

Run:
````sh
docker run --rm   -e USERNAME="$USERNAME"   -e PASSWORD="$PASSWORD"   -e GITHUB_WORKSPACE=/github/workspace   -v "$(pwd)/ena-test:/github/workspace"   emobon-ena-register 2>&1 | tee ena-dev.log
````
