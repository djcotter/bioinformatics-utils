"""
process_bed_file_for_hyena.py
Daniel Cotter

Take in a bed file that contains genomic regions (produced by running pyfaidx --transform bed on a fasta file)
and convert it to the expected input format for training in hyena. Include a flag for --add_metadata which takes a two
column file with names that match the names in the bed file and adds that metadata to the output file.

This script adds the split column to the bed file (split is either train, val, or test) and optionally adds metadata.

When creating the splits, the script ensures that no samples with the same prefix are in different splits (if the --group_by_prefix flag is set).
The prefix is defined as the part of the sample name before the character specified by the --prefix_delimiter flag (default is '_').

The input bed file is expected to have the following columns:
sample_name, start, stop

The output file will have the following columns:
sample_name, start, stop, split, metadata (if --add_metadata is provided)
"""

# imports
import argparse
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Process a bed file for Hyena training."
    )
    parser.add_argument(
        "--input_bed", type=str, required=True, help="Path to the input bed file."
    )
    parser.add_argument(
        "--output_bed", type=str, required=True, help="Path to the output bed file."
    )
    parser.add_argument(
        "--add_metadata",
        type=str,
        default=None,
        help="Path to a two-column metadata file (sample_name, metadata).",
    )
    parser.add_argument(
        "--group_by_prefix",
        action="store_true",
        help="Group samples by prefix when creating splits.",
    )
    parser.add_argument(
        "--prefix_delimiter",
        type=str,
        default="_",
        help="Delimiter to use for extracting prefixes.",
    )
    parser.add_argument(
        "--train_size",
        type=float,
        default=0.7,
        help="Proportion of data to use for training.",
    )
    parser.add_argument(
        "--val_size",
        type=float,
        default=0.15,
        help="Proportion of data to use for validation.",
    )
    parser.add_argument(
        "--test_size",
        type=float,
        default=0.15,
        help="Proportion of data to use for testing.",
    )
    parser.add_argument(
        "--balanced_training",
        action="store_true",
        help="Use balanced training set based on smallest metadata class.",
    )
    return parser.parse_args()


def load_bed_file(bed_path):
    df = pd.read_csv(
        bed_path, sep="\t", header=None, names=["sample_name", "start", "stop"]
    )
    return df


def load_metadata(metadata_path):
    metadata_df = pd.read_csv(
        metadata_path, sep="\t", header=None, names=["sample_name", "metadata"]
    )
    metadata_dict = dict(zip(metadata_df["sample_name"], metadata_df["metadata"]))
    return metadata_dict


def assign_splits(
    df,
    train_size,
    val_size,
    test_size,
    group_by_prefix,
    prefix_delimiter,
    stratify_metadata=False,
    balanced_training=False,
):
    if balanced_training and "metadata" in df.columns:
        # Balanced training set: use train_size of smallest metadata class, then split remainder
        train_list = []
        val_list = []
        test_list = []
        min_class_size = df["metadata"].value_counts().min()
        n_train_per_class = int(train_size * min_class_size)
        val_prop = val_size / (val_size + test_size)
        for meta_class, group in df.groupby("metadata"):
            group = group.sample(frac=1, random_state=42)
            train_samples = group.iloc[:n_train_per_class]
            remaining = group.iloc[n_train_per_class:]
            n_val = int(val_prop * len(remaining))
            val_samples = remaining.iloc[:n_val]
            test_samples = remaining.iloc[n_val:]
            train_samples = train_samples.copy()
            val_samples = val_samples.copy()
            test_samples = test_samples.copy()
            train_samples["split"] = "train"
            val_samples["split"] = "val"
            test_samples["split"] = "test"
            train_list.append(train_samples)
            val_list.append(val_samples)
            test_list.append(test_samples)
        df = pd.concat(train_list + val_list + test_list, ignore_index=True)
        return df
    if group_by_prefix:
        df["prefix"] = df["sample_name"].apply(lambda x: x.split(prefix_delimiter)[0])
        # If stratify_metadata is True and metadata column exists, stratify by metadata class
        if stratify_metadata and "metadata" in df.columns:
            # For each prefix, get its metadata class
            prefix_meta = df.groupby("prefix")["metadata"].first()
            # Stratified split of prefixes by metadata class
            train_prefixes, temp_prefixes = train_test_split(
                prefix_meta.index,
                train_size=train_size,
                stratify=prefix_meta.values,
                random_state=42,
            )
            temp_meta = prefix_meta[temp_prefixes]
            val_prefixes, test_prefixes = train_test_split(
                temp_meta.index,
                test_size=test_size / (test_size + val_size),
                stratify=temp_meta.values,
                random_state=42,
            )
        else:
            unique_prefixes = df["prefix"].unique()
            train_prefixes, temp_prefixes = train_test_split(
                unique_prefixes, train_size=train_size, random_state=42
            )
            val_prefixes, test_prefixes = train_test_split(
                temp_prefixes,
                test_size=test_size / (test_size + val_size),
                random_state=42,
            )

        split_dict = {}
        for prefix in train_prefixes:
            split_dict[prefix] = "train"
        for prefix in val_prefixes:
            split_dict[prefix] = "val"
        for prefix in test_prefixes:
            split_dict[prefix] = "test"

        df["split"] = df["prefix"].map(split_dict)
        df.drop(columns=["prefix"], inplace=True)
    else:
        if stratify_metadata and "metadata" in df.columns:
            train_df, temp_df = train_test_split(
                df, train_size=train_size, stratify=df["metadata"], random_state=42
            )
            val_df, test_df = train_test_split(
                temp_df,
                test_size=test_size / (test_size + val_size),
                stratify=temp_df["metadata"],
                random_state=42,
            )
            train_df["split"] = "train"
            val_df["split"] = "val"
            test_df["split"] = "test"
            df = pd.concat([train_df, val_df, test_df], ignore_index=True)
        else:
            df["split"] = np.random.choice(
                ["train", "val", "test"],
                size=len(df),
                p=[train_size, val_size, test_size],
            )
    return df


def main():
    args = parse_arguments()

    # Load bed file
    bed_df = load_bed_file(args.input_bed)

    # Load metadata if provided
    if args.add_metadata:
        metadata_dict = load_metadata(args.add_metadata)
        bed_df["metadata"] = bed_df["sample_name"].map(metadata_dict)
        # remove rows with NA or empty metadata
        bed_df = bed_df.dropna(subset=["metadata"])
        bed_df = bed_df[bed_df["metadata"] != ""]

    # Assign splits
    bed_df = assign_splits(
        bed_df,
        args.train_size,
        args.val_size,
        args.test_size,
        args.group_by_prefix,
        args.prefix_delimiter,
        stratify_metadata=bool(args.add_metadata),
        balanced_training=args.balanced_training,
    )

    # Save to output file
    # the output order is sample_name, start, stop, split, metadata (if exists)
    if "metadata" in bed_df.columns:
        bed_df = bed_df[["sample_name", "start", "stop", "split", "metadata"]]
    else:
        bed_df = bed_df[["sample_name", "start", "stop", "split"]]
    bed_df.to_csv(args.output_bed, sep="\t", index=False, header=False)
    print(f"Processed bed file saved to {args.output_bed}")


if __name__ == "__main__":
    main()
