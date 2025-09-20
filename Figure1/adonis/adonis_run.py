import argparse
import qiime2 as q2
import pandas as pd
import adonis
import logging
from pathlib import Path

def main(dm_name: str):
    outdir = Path("adonis_out") / dm_name
    outdir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        filename=str(outdir / f"adonis_{dm_name}.log"),
        format='%(asctime)s %(message)s',
        filemode='w',
        level=logging.DEBUG
    )
    logger = logging.getLogger(dm_name)

    dm_paths = {
        'uwuf': 'results/WoL2_3country_uwuf_dm.qza',
        'wuf': 'results/WoL2_3country_wuf_dm.qza',
        'rpca': 'results/rpca_dm.qza',
        'prpca': 'results/WoL2_3country_prpca_dm.qza',
    }
    if dm_name not in dm_paths:
        raise ValueError(f"Unknown dm: {dm_name}")

    dm_art = q2.Artifact.load(dm_paths[dm_name])
    md_q2 = q2.Metadata.load('/home/lakhatib/3country/final_scripts/data/md_for_adonis.tsv')
    vars_to_consider = pd.read_csv(
        '/home/lakhatib/3country/final_scripts/adonis/cols_for_examination.csv',
        index_col=0
    )['columns'].values

    all_res = []
    for v in vars_to_consider:
        if v == 'thdmi_cohort':
            formula = 'thdmi_cohort'
        else:
            formula = f'thdmi_cohort + {v}'

        logger.info(f"[{dm_name}] formula: {formula}")
        res = adonis.adonis_and_reformat(dm_art, md_q2, formula, str(outdir))
        all_res.append(res)

    final = pd.concat(all_res, ignore_index=True)
    final['dm'] = dm_name
    final.to_csv(outdir / f'adonis_results_{dm_name}.csv', index=False)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--dm", required=True)
    args = p.parse_args()
    main(args.dm)
