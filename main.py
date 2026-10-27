"""Run the first project ascent segment and export ideal truth for C."""
import argparse
import json
from pathlib import Path
from ascent import load_config,run_ascent,write_outputs


def main():
    root=Path(__file__).parent
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=root/'config'/'powered_ascent.json')
    parser.add_argument('--output-dir',type=Path,default=None)
    parser.add_argument('--step-s',type=float,default=None)
    args=parser.parse_args()
    config=load_config(args.config)
    name={'early_ascent_attitude_check':'attitude_check','shutdown_separation_check':'shutdown_separation'}.get(config['phase'],'powered_ascent')
    output_dir=args.output_dir or root/'results'/name
    if config['phase']=='shutdown_separation_check':
        from phase_sequence import run_sequence,write_sequence_outputs
        rows,summary=run_sequence(config,step_s=args.step_s)
        write_sequence_outputs(rows,summary,config,output_dir)
    else:
        rows,summary=run_ascent(config,step_s=args.step_s)
        write_outputs(rows,summary,config,output_dir)
    keys=['time_s','total_mass_kg','stage1_propellant_kg','base_altitude_m','airspeed_mps','thrust_n','drag_n']
    print(json.dumps({'phase':summary['phase'],'samples':summary['samples'],
                     'stop_reason':summary['stop_reason'],'final':{k:summary['final'][k] for k in keys}},indent=2))
    print(f'Wrote truth_for_c.csv and source/config metadata to {output_dir}')


if __name__=='__main__':main()
