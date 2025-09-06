#!/usr/bin/env python3
import argparse
import logging
import os
import subprocess
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

def parse_args():
    parser = argparse.ArgumentParser()

    # Archivos de entrada
    parser.add_argument("-p", "--phenotype_file", required=True, help="Path to the (original) phenotype file.")
    parser.add_argument("--phenotypes", default=None, nargs="+")
    parser.add_argument("--relatedness_file") 
    parser.add_argument("--keep_file", default=None, help="Optional file with list of individuals (FID IID) to keep")
    parser.add_argument("--bgen_sample_file") 
    parser.add_argument("--output_dir")
    parser.add_argument("--tmpdir", default="tmp/")
    parser.add_argument("--seed", default=[42], nargs="+", type=int)
    parser.add_argument("--gzip_output", default=False, action="store_true")
    parser.add_argument("--na_code", default="-999")
    parser.add_argument("-o", "--output_file_prefix", required=True)
    parser.add_argument("--overwrite_output", default=False, action="store_true")
    parser.add_argument("--split", action="store_true", help="If set, partition into discovery and replication sets. Otherwise keep all together.")

    return parser.parse_args()


def load_samples(sample_file: str) -> pd.DataFrame:
    df = pd.read_csv(sample_file, sep="\s+", skiprows=2,
                     names=["id_1", "id_2", "missing"])
    df = df.rename(columns={"id_1": "ID"})
    df["ID"] = df["ID"].astype(str)
    return df


def run_greedy_related(rel_df: pd.DataFrame, ids: list, tmpdir: str, frac: float=None, exec_path="GreedyRelated") -> list:
    subset = rel_df.query("ID1 in @ids and ID2 in @ids").drop(
        columns=["HetHet", "IBS0"], errors="ignore"
    )
    
    if frac is None:
        inter_file = f"{tmpdir}/GreedyRelated/greedy_all.dat"
        excl_file = f"{tmpdir}/GreedyRelated/exclude_subjects.txt"
    else:
        part_suffix = f"{int(100*frac)}n{int(100*(1-frac))}"
        inter_file = f"{tmpdir}/GreedyRelated/greedy_{part_suffix}.dat"
        excl_file = f"{tmpdir}/GreedyRelated/exclude_subjects_{part_suffix}.txt"

    os.makedirs(os.path.dirname(inter_file), exist_ok=True)
    subset.to_csv(inter_file, sep="\t", index=False)

    cmd = f"{exec_path} -r {inter_file} --seed 1 --id1 ID1 --id2 ID2 -f Kinship -t 0.0884 -o {excl_file}"
    print("Running:\n\t", cmd)
    subprocess.run(cmd, shell=True, check=True)

    exclude_subjects = pd.read_csv(excl_file, header=None, sep='\t')[0].astype(str).tolist()
    return list(set(ids) - set(exclude_subjects))


def prepare_phenotypes(phenotype_file: str, na_code: str) -> pd.DataFrame:
    df = pd.read_csv(phenotype_file, sep=",")    
    assert "ID" in df.columns, "The phenotype file must contain an 'ID' column. First row is: " + str(df.columns)
    dup_mask = df["ID"].duplicated()
    df.loc[dup_mask, :] = na_code
    df.loc[dup_mask, "ID"] = [-100 * (i + 1) for i in range(dup_mask.sum())]
    df["ID"] = df["ID"].astype(str)
    return df.set_index("ID")


def save_pheno(ids: list, samples_df: pd.DataFrame, pheno_df: pd.DataFrame,
               out_file: str, na_code: str, gzip: bool):
    df = pheno_df.loc[[str(i) for i in ids]].reset_index()
    df = samples_df.merge(df, on="ID", how="left").drop(
        columns=["id_2", "missing"], errors="ignore"
    )
    df.to_csv(out_file, sep="\t", index=False, na_rep=na_code)
    if gzip:
        subprocess.run(f"gzip {out_file}", shell=True)


def main():
    args = parse_args()

    args.bgen_sample_file = os.path.expanduser(args.bgen_sample_file)
    args.output_dir = os.path.expanduser(args.output_dir)

    GreedyExec = "GreedyRelated"

    samples_df = load_samples(args.bgen_sample_file)
    all_ids = samples_df["ID"].tolist()

    rel_df = pd.read_csv(args.relatedness_file, sep="\s+")
    rel_df["ID1"] = rel_df["ID1"].astype(str)
    rel_df["ID2"] = rel_df["ID2"].astype(str)

    pheno_df = prepare_phenotypes(args.phenotype_file, args.na_code)

    os.makedirs(args.tmpdir, exist_ok=True)

    if args.keep_file:
        keep_df = pd.read_csv(args.keep_file, sep="\s+", header=None, names=["ID"])
        keep_ids = set(keep_df["ID"].astype(str))
        samples_df = samples_df[samples_df["ID"].isin(keep_ids)]
        all_ids = samples_df["ID"].tolist()
        logging.info(f"Restricting to {len(all_ids)} individuals from keep file")
 

    if args.split:
        
        fracs_replication = [0.0]
        SEEDS = args.seed

        for SEED in SEEDS:
            out_seed_dir = f"{args.output_dir}/seed_{SEED}/"
            os.makedirs(out_seed_dir, exist_ok=True)
            for frac in fracs_replication:
                np.random.seed(SEED)

                ids_discovery = samples_df.sample(frac=1-frac)["ID"].tolist()
                ids_discovery = run_greedy_related(rel_df, ids_discovery, frac, args.tmpdir, GreedyExec)

                ids_replication = list(set(all_ids) - set(ids_discovery))
                ids_replication = run_greedy_related(rel_df, ids_replication, frac, args.tmpdir, GreedyExec)

                part_suffix = f"{int(100*frac)}n{int(100*(1-frac))}"

                pd.DataFrame(ids_discovery).to_csv(
                    f"{args.tmpdir}/GreedyRelated/ids_discovery_{part_suffix}.txt",
                    index=False, header=False)
                pd.DataFrame(ids_replication).to_csv(
                    f"{args.tmpdir}/GreedyRelated/ids_replication_{part_suffix}.txt",
                    index=False, header=False)

                disc_file = f"{out_seed_dir}/{args.output_file_prefix}-discovery_{part_suffix}.csv"
                repl_file = f"{out_seed_dir}/{args.output_file_prefix}-replication_{part_suffix}.csv"

                save_pheno(ids_discovery, samples_df, pheno_df, disc_file, args.na_code, args.gzip_output)
                save_pheno(ids_replication, samples_df, pheno_df, repl_file, args.na_code, args.gzip_output)

                logging.info(f"Seed={SEED}, frac={frac}: {len(ids_discovery)} discovery, {len(ids_replication)} replication")
    else:

        out_dir = args.output_dir
        ids_final = run_greedy_related(rel_df, all_ids, tmpdir=args.tmpdir, exec_path=GreedyExec)

        out_file = f"{out_dir}/{args.output_file_prefix}.csv"
        save_pheno(ids_final, samples_df, pheno_df, out_file, args.na_code, args.gzip_output)

        pd.DataFrame(ids_final).to_csv(f"{args.tmpdir}/GreedyRelated/ids_all.txt", index=False, header=False)

        logging.info(f"Seed={args.seed}, no split: {len(ids_final)} individuals retained after relatedness filtering")
        logging.info(f"Phenotype file saved to {out_file=}")

    # for SEED in SEEDS:

    #     out_seed_dir = f"{args.output_dir}/seed_{SEED}/"
    #     os.makedirs(out_seed_dir, exist_ok=True)

    #     for frac in fracs_replication:

    #         np.random.seed(SEED)

    #         ids_discovery = samples_df.sample(frac=1-frac)["ID"].tolist()
    #         ids_discovery = run_greedy_related(rel_df, ids_discovery, frac, args.tmpdir, GreedyExec)

    #         ids_replication = list(set(all_ids) - set(ids_discovery))
    #         ids_replication = run_greedy_related(rel_df, ids_replication, frac, args.tmpdir, GreedyExec)

    #         part_suffix = f"{int(100*frac)}n{int(100*(1-frac))}"

    #         pd.DataFrame(ids_discovery).to_csv(
    #             f"{args.tmpdir}/GreedyRelated/ids_discovery_{part_suffix}.txt",
    #             index=False, header=False)
    #         pd.DataFrame(ids_replication).to_csv(
    #             f"{args.tmpdir}/GreedyRelated/ids_replication_{part_suffix}.txt",
    #             index=False, header=False)

    #         disc_file = f"{out_seed_dir}/{args.output_file_prefix}-discovery_{part_suffix}.csv"
    #         repl_file = f"{out_seed_dir}/{args.output_file_prefix}-replication_{part_suffix}.csv"

    #         save_pheno(ids_discovery, samples_df, pheno_df, disc_file, args.na_code, args.gzip_output)
    #         save_pheno(ids_replication, samples_df, pheno_df, repl_file, args.na_code, args.gzip_output)

    #         logging.info(f"Seed={SEED}, frac={frac}: {len(ids_discovery)} discovery, {len(ids_replication)} replication")


if __name__ == "__main__":
    main()
