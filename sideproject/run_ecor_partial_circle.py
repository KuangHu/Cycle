"""Run partial circle detection on all ECOR samples × references."""
import os, sys, logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
sys.path.insert(0, '/home/kuangh/tools/Cycle')

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(name)s: %(message)s')
logger = logging.getLogger(__name__)


def run_one(ref, sample_id, guide_json, fastq, pc_dir):
    import logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(name)s: %(message)s')
    from cycle.partial_circle.detector import PartialCircleDetector
    from pathlib import Path
    detector = PartialCircleDetector(output_dir=str(pc_dir), threads=2)
    detector.run_sample(
        is_records_path=Path(guide_json),
        sample_id=sample_id,
        fastq_path=Path(fastq),
    )
    return f"{ref}/{sample_id}"


BASE = Path("/groups/rubin/projects/kuang/out/IS110/ECOR_batch/output")
FASTQ_DIR = Path("/groups/rubin/projects/kuang/out/IS110/ECOR_batch/fastq")
fastq_map = {fq.name.split(".")[0]: fq for fq in FASTQ_DIR.glob("*.fastq.gz")}

tasks = []
for ref in ["K12_MG1655", "CFT073", "IAI39", "Sakai"]:
    ref_dir = BASE / ref
    if not ref_dir.is_dir():
        continue
    for sample_dir in sorted(ref_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        sample_id = sample_dir.name
        guide_json = sample_dir / "is_formatter_output" / sample_id / f"{sample_id}_is_records_guide.json"
        if not guide_json.exists():
            continue
        fastq = fastq_map.get(sample_id)
        if not fastq:
            continue
        pc_dir = sample_dir / "partial_circle_output"
        pc_summary = pc_dir / sample_id / f"{sample_id}_partial_circle_summary.json"
        if pc_summary.exists():
            logger.info("Skipping %s/%s (done)", ref, sample_id)
            continue
        tasks.append((ref, sample_id, str(guide_json), str(fastq), str(pc_dir)))

logger.info("Running %d tasks with 4 parallel workers", len(tasks))

with ProcessPoolExecutor(max_workers=4) as pool:
    futures = {pool.submit(run_one, *t): t[1] for t in tasks}
    done = 0
    for fut in as_completed(futures):
        done += 1
        try:
            result = fut.result()
            logger.info("  %d/%d done: %s", done, len(tasks), result)
        except Exception as e:
            logger.error("  %d/%d FAILED %s: %s", done, len(tasks), futures[fut], e)

logger.info("=== Done ===")
