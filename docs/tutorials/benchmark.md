# HADDOCK3 benchmark tutorial

Developing HADDOCK3 requires handling thousands of parameters, either at
the CNS level or the Python level. Therefore, unit-testing is sometimes
not enough. To assess the correct functioning of HADDOCK3 we have
assembled a benchmarking workflow where we dock well known and studied
complexes and ensure the outcome matches expectations.

To run HADDOCK3 as a user, you don't need to run this benchmark. But if
you are an advanced user who has configured HADDOCK3 in your systems,
you may wish to run the benchmark to make sure HADDOCK3 behaves as
expected.

HADDOCK3 has two command-line clients that set up the benchmark workflow
automatically. Let's go through them step-by-step.

## Obtaining the initial models

The Bonvin lab has manually curated a list of systems ideal for this
benchmark. The repository is available here. Clone it to your machine
with:

`git clone https://github.com/haddocking/BM5-clean`

If you inspect the new `BM5-clean` folder, you will see the
`HADDOCK-ready` folder has more than 200 sub-folders with "ready-to-use"
systems for HADDOCK.

## Setting up the runs

Having HADDOCK3 installed, you can set up the runs with the following
command:

```
# first, make a new folder where you want the runs to be stored
mkdir <folder-where-the-BM-jobs-will-be-created>
haddock3-bm BM5-clean/HADDOCK-ready <folder-where-the-BM-jobs-will-be-created> [OPTIONS]
```

`[OPTIONS]` are optional :relaxed:. Look at the `haddock3-bm -h` help
menu and select which ones better adapt to your system.

## Running with the daemon

Running more than 200 HADDOCK jobs is time-consuming. Inspecting the
queue system manually over several days is a pain. Therefore, we created
a daemon that will handle the queue and submit the jobs for you. Use the
command:

`haddock3-dmn <folder-where-the-BM-jobs-were-created> [OPTIONS]`

Use `haddock3-dmn -h` to read all possible options.

From this point on, the daemon will manage the queue for you. You can
also run the daemon as a job; please read further.

Initially, all jobs have a file named `AVAILABLE` pointing that the jobs
are ready to run. Jobs that are sent to the queue get the file flag
`RUNNING`. Completed jobs get the file flag `DONE`, and failed jobs get
the file flag `FAIL`.You can inspect which jobs are at which state with
the `find` command:

`find <folder-where-the-BM-jobs-were-created> -name DONE`

You can count the number of folders with:

`find <folder-where-the-BM-jobs-were-created> -name DONE | wc -l`

## Stopping the daemon

If you want to stop the daemon, press `Ctrl+c` in the window where you
launch it. And manually cancel the remaining jobs. You can try the
following command:

`qstat -a | awk '{print $4}' | grep BM5 | xargs scancel`

`BM5` comes because all jobs prepared for the benchmark have the `BM5`
tag in their job name (unless you defined otherwise in the `[OPTIONS]`)

## Further notes

As with any command-line tool, the haddock3-daemon will run as long as
the terminal window is open. However, we don't recommend using the `&`
terminal detacher suffix command unless you know what you do.
Alternatively, you could use [`tmux`](https://github.com/tmux/tmux) to
maintain the session running even when you disconnect from your login
node.

A better way can be to send the daemon job to some queue. You can easily
do that by using the daemon job created by the `haddock3-bm` command.

## What happens if you halt the daemon?

If you halt the daemon, you can restart the operations from where they
left. Two scenarios are possible:

1. you also canceled running jobs
2. you let the running jobs finish properly

In the first scenario, you want the previously running jobs to restart.
For that, run the command:

`haddock3-dmn <new-folder-of-my-preference> --restart [OPTIONS]`

Or you can update the `hd3-daemon.job` file with the `--restart` option.

In the second scenario, run the daemon as you would from scratch. Only
the jobs with the file flag `AVAILABLE` will be executed.

Thanks for using HADDOCK3. [Any feedback is very
welcomed](https://github.com/haddocking/haddock3/issues).

The HADDOCK team.
