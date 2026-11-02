import csv
import hashlib
import json
import math
from pathlib import Path
from urllib.request import urlopen


REVISION = "b245d3b81aa36b7941ec10f3f4b508999d106a6d"
SOURCE = f"https://raw.githubusercontent.com/shahar603/Telemetry-Data/{REVISION}/SpaceX%20CRS-10/JSON/SpaceX%20CRS-10%20raw.json"


def main():
    raw = urlopen(SOURCE, timeout=60).read()
    data = json.loads(raw)
    lengths = {len(data[key]) for key in ("time", "velocity", "altitude")}
    if len(lengths) != 1:
        raise ValueError("Source column lengths differ.")
    samples = [(float(t), float(v), float(h) * 1000)
               for t, v, h in zip(data["time"], data["velocity"], data["altitude"])
               if all(math.isfinite(x) for x in (t, v, h)) and 0 <= t <= 140 and v >= 0 and h >= 0]
    rows = []
    for target in range(141):
        time_s, speed_mps, altitude_m = min(samples, key=lambda row: abs(row[0] - target))
        if abs(time_s - target) > 0.15:
            continue
        rows.append({"target_time_s": target, "time_s": time_s,
                     "speed_mps": speed_mps, "altitude_m": altitude_m})
    output = Path(__file__).resolve().parents[1] / "data" / "crs10"
    output.mkdir(parents=True, exist_ok=True)
    with (output / "ascent_reference.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "classification": "E",
        "origin": "Reprocessed author-published webcast digitization, not a new direct video extraction or onboard telemetry release.",
        "source_url": SOURCE,
        "source_revision": REVISION,
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "reference_video_url": "https://www.youtube.com/watch?v=rUDLxFUMC9c",
        "source_readme": f"https://github.com/shahar603/Telemetry-Data/blob/{REVISION}/README.md",
        "source_mission_readme": f"https://github.com/shahar603/Telemetry-Data/blob/{REVISION}/SpaceX%20CRS-10/README.md",
        "source_columns": {"time": "s", "velocity": "m/s", "altitude": "km"},
        "processing": "Nearest available sample within 0.15 s of each integer target second; actual source time retained; altitude multiplied by 1000; no interpolation or smoothing.",
        "source_rows": len(data["time"]),
        "exported_rows": len(rows),
        "missing_target_seconds": sorted(set(range(141)) - {row["target_time_s"] for row in rows}),
        "coverage": "Only 0..140 s, conservatively before the author's estimated MECO at 143 s; combined stack ascent. No post-separation booster reference.",
        "speed_reference_frame": "Not independently established; not automatically air-relative speed.",
        "limitations": "No new video-frame verification performed. Display quantization, OCR errors and stage mixing remain source risks. Published analysed trajectory/acceleration are not treated as measurements."
    }
    (output / "ascent_reference_source.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"rows": len(rows), "first": rows[0], "last": rows[-1]}, indent=2))


if __name__ == "__main__":
    main()
