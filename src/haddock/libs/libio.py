"""General Input/Output routines."""
import shutil


def save_stream_to_file(stream, file_name):
    stream.seek(0)
    with open(file_name, 'w') as fout:
        shutil.copyfileobj(stream, fout)


def save_streams_to_files(streams, file_names):
    for stream, file_name in zip(streams, file_names):
        save_stream_to_file(stream, file_name)
