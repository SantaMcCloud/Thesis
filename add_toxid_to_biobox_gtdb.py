import subprocess
import argparse
import os


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="",
        description="Use this scrip to add toxid from GTBD to binning files",
        usage="python add_toxid_to_biobox_gtdb.py -b BINNING FILE -g GTDB FILE PATH",
        add_help=False
    )

    parser.add_argument("--gtdb_path", "-g", type=str)
    parser.add_argument("--binning_file", "-b", type=str)
    
    return parser.parse_args()

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

def get_gtdb_names(gtdb_file_path):
    print("Start mapping the classified names to the sequences.")
    mapping = {}
    if not gtdb_file_path.endswith("/"):
        gtdb_file_path+= "/"
    for dir in os.listdir(gtdb_file_path):
        subpath = gtdb_file_path + dir
        for subdir in os.listdir(subpath):
            path = subpath + "/" + subdir
            for file in os.listdir(path):
                print(f"Load file: {dir}/{subdir}/{file}")
                with open("{0}{1}/{2}/{3}".format(gtdb_file_path, dir, subdir, file),"r") as f:
                    for lines in f:
                        if lines.startswith("user_genome"):
                            continue
                        line = lines.split("\t")
                        seqid = line[0]
                        if ";" in line[1]:
                            lineage = line[1].split(";")
                            name = lineage[-1]
                        else:
                            name = line[1]
                        if name not in mapping.keys():
                            mapping[name] = [seqid]
                        else:
                            mapping[name].append(seqid)
    print(f'# Unclassified Archaea: {len(mapping['Unclassified Archaea'])}\n # Unclassified Bacteria: {len(mapping['Unclassified Bacteria'])}\n # Unclassified: {len(mapping['Unclassified'])}')
    return mapping

def generate_ncbi_names(gtdb_names):
    print("Generate the NCBI names with gtdb_to_taxdum.")
    with open("gtdb_names.txt", "w") as file:
        for name in gtdb_names.keys():
            file.write(f'{name}\n')
    arg = ["ncbi-gtdb_map.py", '-q', 'gtdb_taxonomy', 'gtdb_names.txt', 'bac120_metadata_r220.tsv', 'ar53_metadata_r220.tsv']
    subprocess.run(arg)

def generate_taxids():
    names = {}
    print("Extracted NCBI names from ncbi-gtdb/taxonomy_map_summary.tsv so taxonkit then can generate the taxids.")
    with open("ncbi-gtdb/taxonomy_map_summary.tsv", "r") as gtdb, open("ncbi_names.txt", "w") as ncbi:
        for lines in gtdb:
            if 'taxonomy' in lines:
                continue
            line = lines.split('\t')
            ncbi.write(f'{line[1].split("_")[-1]}\n')
            if line[0] not in names.keys():
                names[line[0]] = line[1].split("_")[-1]
    arg = ["taxonkit", "name2taxid", 'ncbi_names.txt']
    with open("taxids.txt", 'w') as file:
        subprocess.run(arg, stdout=file)
    return names

def create_file(binning, gtdb_names, name_map, binning_file):
    taxid = {}
    file_name = binning_file.split(".")
    with open("{0}_GTDB.{1}".format(file_name[0], file_name[1]), "w") as binning_f:
        binning_f.write("#CAMI Format for Binning\n")
        binning_f.write("@Version:0.9.0\n")
        binning_f.write("@SampleID:_SAMPLEID_\n")
        binning_f.write("@@SEQUENCEID\tBINID\tTAXID\n")
        with open("taxids.txt", "r") as tax:
            for lines in tax:
                line = lines.split("\t")
                taxid[line[0]] = line[1].replace("\n", "")
        for gtdb_name in gtdb_names.keys():
            for seqid in gtdb_names[gtdb_name]:
                seqid = seqid.replace("_", "|")
                if seqid in binning.keys():
                    binid = binning[seqid][0]
                    ncbi_name = name_map[gtdb_name]
                    tax = taxid[ncbi_name]
                    binning_f.write(f"{seqid}\t{binid}\t{tax}\n")

def run():
    args = parse_arguments()
    if not args.gtdb_path:
        print("Please use followed command with all arguments!")
        print("python add_toxid_to_biobox_gtdb.py -b BINNING FILE -g GTDB FILE PATH")
    else:
        binning = load_binning(args.binning_file)
        gtdb_names = get_gtdb_names(args.gtdb_path)
        generate_ncbi_names(gtdb_names)
        name_map = generate_taxids()
        create_file(binning, gtdb_names, name_map, args.binning_file)

if __name__=="__main__":
    run()