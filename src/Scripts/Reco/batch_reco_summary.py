import json, os, argparse
import dask
import awkward as ak
from coffea import processor
from coffea.nanoevents import BaseSchema
from coffea.dataset_tools import (
    apply_to_fileset,
    max_chunks,
    preprocess,
)
import uproot
import logging

# Configure logger
logging.basicConfig(level=logging.INFO, format="%(asctime)s — %(levelname)s — %(message)s")
logger = logging.getLogger(__name__)
for noisy in ["uproot", "dask", "fsspec"]:
    logging.getLogger(noisy).setLevel(logging.WARNING)


def remove_empty_files(fileset):
    cleaned = {}
    for dataset, content in fileset.items():
        files = content.get('files', {})
        valid = {}
        for fp, tn in files.items():
            try:
                with uproot.open(f"{fp}:{tn}") as tree:
                    if tree.num_entries > 0:
                        valid[fp] = tn
            except Exception as e:
                logger.warning(f"Skipping bad file {fp}: {e}")
        if valid:
            cleaned[dataset] = {'files': valid}
    return cleaned


class CountProcessor(processor.ProcessorABC):
    def process(self, events):
        total = int(ak.num(events, axis=0).compute())
        bad = int(ak.sum(events.chi2_status != 0).compute())
        return {"entries": total, "bad_reco": bad}

    def postprocess(self, accumulator):
        return accumulator


def main():
    parser = argparse.ArgumentParser(description="Batch reco summary for all MC processes")
    parser.add_argument('-e','--era', required=True)
    parser.add_argument('-t','--tag', required=True)
    args = parser.parse_args()

    json_path = f'/home/mukund/Projects/PhysicsTools/NanoAODTools/Datasets/{args.tag}_reco_{args.era}_dataFiles.json'
    with open(json_path, 'r') as f:
        dicti = json.load(f)

    fileset_all = {}
    for proc_name, files in dicti.get('MC_mu', {}).items():
        datasetName = f"{args.era}_{proc_name}"
        fileset_all[datasetName] = { 'files': files }

    fileset_all = remove_empty_files(fileset_all)
    logger.info(f"Found {len(fileset_all)} cleaned datasets")

    if not fileset_all:
        logger.error("No datasets to run")
        return

    dataset_runnable, dataset_updated = preprocess(fileset_all, align_clusters=False, files_per_batch=1, skip_bad_files=True, save_form=False)

    to_compute = apply_to_fileset(CountProcessor(), max_chunks(dataset_runnable, 300), schemaclass=BaseSchema)
    (out,) = dask.compute(to_compute, scheduler='threads')

    # out is dict: dataset -> accumulator
    print("Dataset summary:\n")
    for d, acc in out.items():
        entries = acc.get('entries', 0)
        bad = acc.get('bad_reco', 0)
        print(f"{d}: entries={entries}, bad_reco={bad}")

    # Optionally save summary
    out_path = f"outputs/{args.era}_batch_reco_summary_{args.tag}.json"
    os.makedirs('outputs', exist_ok=True)
    with open(out_path, 'w') as of:
        json.dump(out, of, default=int, indent=2)
    logger.info(f"Saved summary to {out_path}")


if __name__ == '__main__':
    main()
