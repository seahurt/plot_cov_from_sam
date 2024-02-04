#!/usr/bin/env python3

import logging.handlers
import os
from pathlib import Path

import fire
import matplotlib.ticker as mtick
import pandas
import seaborn as sns
from matplotlib import pyplot as plt

logger = logging.getLogger(__name__)


def plot_sam(sam, mode, outfile=None, title=None, _logger=logger, dpi=100,
             bin_count=500):
    sam = str(sam)
    if not title:
        title = Path(sam).name
    if not outfile:
        outfile = sam + ".png"
    if mode == "coverage":
        plot_hist_sam(sam, title, outfile, _logger=_logger, dpi=dpi,
                      bin_count=bin_count)
    elif mode == "coverage3":
        plot_line_sam(sam, title, outfile, _logger=_logger, dpi=dpi)
    else:
        raise ValueError(f"unsupported mode: {mode}")


def plot_hist_sam(samfile, title, outfile, _logger=logger, dpi=100,
                  bin_count=500):
    """get total length and seq id"""

    acc2length, total_length = parse_sam_sq(samfile)

    offset = 0
    acc2offset = {}
    for acc, length in acc2length.items():
        acc2offset[acc] = offset
        offset += length
    df = make_bin_data(samfile, acc2offset, total_length, bin_count)
    df.to_csv(f"{samfile}.csv")
    plot_cov(df, title, outfile, bin_count=bin_count, _logger=_logger, dpi=dpi)


def plot_cov(df, title, outfile, bin_count, _logger=logger, dpi=100):
    _logger.info(f"plot {outfile}")
    if df is None:
        _logger.warning(f"df is none, won't plot")
        return

    if len(df) == 0:
        _logger.warning("no reads aligned, won't plot")
        return

    sns.set_theme(style="darkgrid")
    plt.subplots(figsize=(10, 5))
    sns.histplot(list(df['bin_index']), stat='count', binwidth=1, binrange=(0, bin_count))
    plt.xlabel('Coverage (%)', fontweight='bold')
    plt.ylabel('Reads Count', fontweight='bold')
    # plt.xlim(0, bin_count+1)
    plt.xlim(0, bin_count+1)
    # 设置X主坐标label显示百分比
    plt.gca().xaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: str(int(x / bin_count * 100))))

    plt.title(title, fontweight='bold')
    plt.tight_layout()
    plt.savefig(outfile, format='png', dpi=dpi)
    plt.close()


def plot_line_sam(samfile, title, outfile, _logger=logger, dpi=100):
    acc2length, total_length = parse_sam_sq(samfile)

    offset = 0
    acc2offset = {}
    for acc, length in acc2length.items():
        acc2offset[acc] = offset
        offset += length

    df = make_depth_data(samfile, acc2offset, total_length)
    df.to_csv(f"{samfile}.csv")
    plot_cov2(df, title, outfile, _logger=_logger, dpi=dpi)


def plot_cov2(df, title, outfile, _logger=logger, dpi=100):
    _logger.info(f"plot {outfile}")
    if df is None:
        _logger.warning(f"df is none, won't plot")
        return

    if len(df) == 0:
        _logger.warning("no reads aligned, won't plot")
        return
    sns.set_theme(style="darkgrid")
    plt.subplots(figsize=(10, 5))
    sns.lineplot(df, x="position_ratio", y="depth")
    plt.xlabel('Coverage(%)', fontweight='bold')
    plt.ylabel('Depth', fontweight='bold')
    plt.xlim(0, 101)
    plt.ylim(0)
    # plt.xticks(np.arange(0, bin_count+1, step=0.2*bin_count))
    plt.title(title, fontweight='bold')
    plt.tight_layout()
    plt.savefig(outfile, format='png', dpi=dpi)
    plt.close()


def parse_sam_sq(sam):
    acc2length = {}
    with open(sam) as f:
        for line in f:
            if line.startswith('@SQ'):
                parts = line.strip().split()
                acc = parts[1].split(':', 1)[1]
                length = int(parts[2].split(':')[1])
                acc2length[acc] = length
            elif not line.startswith('@'):
                break
    return acc2length, sum(acc2length.values())


def get_read_length(samfile):
    with open(samfile) as f:
        for line in f:
            if line.startswith('@'):
                continue
            info = line.split()
            if len(info) < 12:
                continue
            return len(info[9])


def make_bin_data(samfile, acc2offset: dict, total_ref_length, bin_count):
    step = total_ref_length // bin_count
    if step == 0:
        step = 1

    records = []
    with open(samfile) as f:
        for line in f:
            if line.startswith('@'):
                continue
            info = line.split()
            acc = info[2]
            if acc == "*":
                continue
            pos = int(info[3])
            abs_pos = pos + acc2offset[acc]
            records.append({'accession': acc,
                            'position': pos,
                            'abs_position': abs_pos,
                            'bin_index': abs_pos // step
                            })
    records.sort(key=lambda x: x['abs_position'])

    df = pandas.DataFrame.from_records(records)
    return df


def make_depth_data(samfile, acc2offset, total_ref_length):
    samfile = Path(samfile)
    cmd = f"samtools sort {samfile} -o {samfile}.bam && samtools depth {samfile}.bam -aa > {samfile}.depth && rm -rf {samfile}.bam"
    ret = os.system(cmd)
    if ret != 0:
        raise ValueError(f"calc depth err for {samfile}")
    df = pandas.read_csv(str(samfile) + ".depth", header=None, sep="\t")
    df.columns = ["accession", "position", "depth"]
    df["offset"] = [acc2offset.get(acc, 0) for acc in df["accession"]]
    df["abs_position"] = df["position"] + df["offset"]
    df["position_ratio"] = df["abs_position"] / total_ref_length * 100
    return df


if __name__ == '__main__':
    fire.Fire(plot_sam)
