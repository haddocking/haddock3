# Benchmark Deamon
#
#    (_) L|J
#    (")  |
#    /_\--|
#  _/\ /  |
#    _W_  |
#
import argparse
import glob
import subprocess
import os
import time


class Job:

    def __init__(self, job_f):
        self.job_filename = job_f
        self.job_path = '/' + '/'.join(self.job_filename.split('/')[1:-1])
        self.done_check = f'{self.job_path}/DONE'
        self.running_check = f'{self.job_path}/RUNNING'
        self.fail_check = f'{self.job_path}/FAIL'
        self.status = None

    def get_status(self):

        if os.path.isfile(self.done_check):
            self.status = 'DONE'
        elif os.path.isfile(self.running_check):
            self.status = 'RUNNING'
        elif os.path.isfile(self.fail_check):
            self.status = 'FAIL'
        else:
            self.status = 'AVAILABLE'

        return self.status

    def submit(self):
        subprocess.run(['qsub', self.job_filename])


class Queue:

    def __init__(self, job_list, job_limit):
        self.job_list = [Job(j) for j in job_list]
        self.job_limit = job_limit

    def get_empty_spots(self):
        return self.job_limit - self.current_jobs()

    def list_available_jobs(self):
        avail = []
        for j in self.job_list:
            status = j.get_status()
            if status == 'AVAILABLE':
                avail.append(j)
        return avail

    @staticmethod
    def current_jobs():
        concurrent_cmd = "qstat -a | awk '{print $4}' | grep BM5 | wc -l"
        p = subprocess.Popen(concurrent_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = p.communicate()[0]
        return int(output.decode('utf-8'))


if __name__ == '__main__':
    ''' Deamon to submit/check all Benchmark jobs '''

    parser = argparse.ArgumentParser(description='This is the submission Daemon')
    parser.add_argument("benchmark_path",
                        help="Location of BM5")
    parser.add_argument("--job_limit",
                        help="How many jobs should run at the same time",
                        action="store",
                        default=10)

    args = parser.parse_args()

    path = args.benchmark_path
    job_limit = int(args.job_limit)
    job_list = glob.glob(f'{path}/*/*.job')
    if not job_list:
        print('+ ERROR: No jobs found')
        exit()
    job_list.sort()

    queue = Queue(job_list, job_limit)

    while 1:

        # submit until the empty spots are over
        available_jobs = queue.list_available_jobs()
        empty_spots = queue.get_empty_spots()

        if available_jobs and empty_spots != 0:
            for job in available_jobs:
                if empty_spots > 0:
                    job.submit()
                    time.sleep(5)
                    empty_spots -= 1
                else:
                    break
        elif not available_jobs:
            # if there are no available jobs, there is nothing else to be done!
            exit()

        # chill
        time.sleep(120)
