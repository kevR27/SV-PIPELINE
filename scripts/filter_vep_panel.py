#!/usr/bin/env python3

import argparse


def load_genes(path):
    with open(path) as f:
        return {
            line.strip()
            for line in f
            if line.strip() and not line.startswith("#")
        }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vep", required=True)
    parser.add_argument("--genes", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    genes = load_genes(args.genes)

    with open(args.vep) as fin, open(args.output, "w") as fout:
        csq_fields = None

        for line in fin:

            if line.startswith("##INFO=<ID=CSQ"):
                fout.write(line)

                if "Format:" in line:
                    fmt = line.split("Format:", 1)[1]
                    fmt = fmt.split('">', 1)[0]
                    csq_fields = fmt.split("|")

                continue

            if line.startswith("#"):
                fout.write(line)
                continue

            if csq_fields is None:
                continue

            symbol_idx = (
                csq_fields.index("SYMBOL")
                if "SYMBOL" in csq_fields
                else None
            )

            if symbol_idx is None:
                continue

            info = line.rstrip("\n").split("\t")[7]

            matched = False

            for field in info.split(";"):
                if not field.startswith("CSQ="):
                    continue

                csq = field[4:]

                for annotation in csq.split(","):
                    values = annotation.split("|")

                    if symbol_idx < len(values):
                        symbol = values[symbol_idx]

                        if symbol in genes:
                            matched = True
                            break

                if matched:
                    break

            if matched:
                fout.write(line)


if __name__ == "__main__":
    main()
