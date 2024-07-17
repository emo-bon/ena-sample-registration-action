The action that will register sample accessions in ENA. Input to this comes from 
* The files in the [batches](https://github.com/emo-bon/sequencing-data/tree/main/shipment) folders
* The [checklist translations](https://github.com/emo-bon/sequencing-profile)

The action creates XML files which are then uploaded to ENA, returning sample accession numbers, which are then put back into the relevant files in the [batches](https://github.com/emo-bon/sequencing-data/tree/main/shipment) folders.
The code that sits behind this is called action.py
