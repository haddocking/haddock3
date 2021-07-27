from pathlib import Path
import re
import sys
import argparse
import logging

logging.basicConfig(level="DEBUG",
                    format=("[%(asctime)s] L%(lineno)d %(levelname)s - %(message)s"),
                    datefmt="%d/%m/%Y %H:%M:%S")


def retrieve_dependency(cns_script):
    """Read a CNS script and identify its dependencies."""
    dep_list = []
    marker_regex = r"@RUN:(.*\.cns|.*\.inp)"
    with open(cns_script, 'r') as cns_fh:
        for line in cns_fh.readlines():
            match = re.findall(marker_regex, line)
            if match:
                dependency = match[0].split('/')[-1]
                dep_list.append(dependency)
    return list(set(dep_list))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('cns_script')
    parser.add_argument('protocols_folder')
    args = parser.parse_args()

    protocol_dir = Path(args.protocols_folder)
    if not protocol_dir.exists():
        logging.error(f"{protocol_dir} does not exist.")
        sys.exit()

    main_cns_script = Path(args.cns_script)
    if not main_cns_script.exists():
        logging.error(f"{main_cns_script} does not exist")
        sys.exit()

    logging.info(f'Protocols directory: {protocol_dir}')
    logging.info(f'Retrieving dependencies of {main_cns_script}')

    dependency_list = retrieve_dependency(main_cns_script)

    if dependency_list:
        logging.info(f'First level dependencies: {len(dependency_list)}')

    logging.info('Looking for deeper dependencies')
    done = False
    evaluated_dic = dict((k, False) for k in dependency_list)
    while not done:
        for cns_script in dependency_list:
            if cns_script in evaluated_dic:
                if evaluated_dic[cns_script]:
                    continue

            loc = protocol_dir / cns_script
            if not loc.exists():
                logging.error(f'Dependency not found, investigate: {loc}')
                sys.exit()

            dep_list = retrieve_dependency(loc)
            for dependency in dep_list:
                if dependency not in evaluated_dic:
                    evaluated_dic[dependency] = False
                    dependency_list.append(dependency)

            evaluated_dic[cns_script] = True

        if all(evaluated_dic.values()):
            done = True

    # make a copy friendly string
    total_dependencies = len(evaluated_dic.keys())
    logging.info(f'Total dependencies found: {total_dependencies}')

    dependency_str = ' '.join(map(str, evaluated_dic.keys()))
    logging.info(f'{main_cns_script.name} depends on: {dependency_str}')
