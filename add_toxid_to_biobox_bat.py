import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="",
        description="Use this scrip to add the toxid from BAT to the binning file",
        usage="python add_toxid_to_biobox_bat.py -b BINNING FILE -p ORF FILE",
        add_help=False
    )

    parser.add_argument("--binning_file", "-b", type=str)
    parser.add_argument("--orf_file", "-p", type=str)
    parser.add_argument("--semi", "-s", action='store_true')

    args = parser.parse_args()

    return args

def load_binning(binning_file):
    print("Start with extracting of the binning file.")
    print(f"Load {binning_file}")
    binning = {}
    with open(binning_file, "r") as file:
        for line in file:
            if line.startswith("@") or line.startswith("#"):
                continue
            line = line.replace("\n", "")
            seqid, binid = line.split("\t")
            binning[seqid] = (binid, None)
    return binning

def load_orf_files(orf_file, semi):
    print("Start with extracting the taxid from the ORF file.")
    print(f"Load {orf_file}")
    orf = {}
    with open(orf_file, "r") as bin:
        for line in bin:
            if line.startswith("#"):
                continue
            line_split = line.split("\t")
            if semi:
                id = line_split[0].split("_")[2]
            else:
                id = line_split[0].split("_")[1]
            lineage = line_split[3].split(";")
            if len(lineage) == 1 and lineage[0] == "":
                continue
            if id not in orf:
                orf[id] = lineage
            else:
                if len(orf[id]) < len(lineage):
                    orf[id] = lineage
    return orf

def merge_dict(binning, orf):
    print("Merge the dict to be able to write the binning file with toxid!")
    for key in binning.keys():
        toxid = orf[key][-1]
        if "*" in toxid:
            toxid = toxid[:-1]
        binning[key] = (binning[key][0], toxid)
    return binning

def create_file(final_dict, binning_file):
    print("Create binning file with TOXID column")
    name = binning_file.split(".")
    with open("{0}_BAT.{1}".format(name[0], name[1]), "w") as binning:
        binning.write("#CAMI Format for Binning\n")
        binning.write("@Version:0.9.0\n")
        binning.write("@SampleID:_SAMPLEID_\n")
        binning.write("@@SEQUENCEID\tBINID\tTAXID\n")
        for key in final_dict.keys():
            id = key
            binid, taxid = final_dict[key]
            binning.write(f"{id}\t{binid}\t{taxid}\n")

def run():
    args = parse_arguments()
    if not args.binning_file or not args.orf_file:
        print("Please use followed command with all arguments!")
        print("python add_toxid_to_biobox_bat.py -b BINNING FILE -p ORF FILE")
    else:
        binning = load_binning(args.binning_file)
        orf = load_orf_files(args.orf_file, args.semi)
        final_dict = merge_dict(binning, orf)
        create_file(final_dict, args.binning_file)

if __name__ == "__main__":
    run()