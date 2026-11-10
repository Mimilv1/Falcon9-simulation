from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys
import textwrap
import zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / 'report'
sys.path.insert(0, str(ROOT))
from frames import body_to_ecef_matrix
from VehicleData import build_falcon9

BLUE, GREEN, ORANGE, GRAY = '#125B8E', '#237567', '#AB6422', '#7D8790'
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':8,
    'axes.spines.top':False, 'axes.spines.right':False, 'pdf.fonttype':42})
FILES = ['truth.csv', 'stage2_truth.csv', 'navigation.csv', 'measurements.csv',
    'events.json', 'stage2_events.json', 'event_timeline.csv', 'validation.csv',
    'summary.json', 'navigation_covariance.npz', 'startracker_innovations.csv']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_inputs(run, canonical):
    summary=json.loads((run/'summary.json').read_text())
    normalized={}; changed_endings=[]
    for name,sha in summary['input_sha256'].items():
        source=ROOT/name
        normalized[name]=hashlib.sha256(source.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
        if digest(source)!=sha:
            if canonical.get(name)!=normalized[name]:
                raise ValueError('Runtime source differs: '+name)
            changed_endings.append(name)
    return normalized,changed_endings


def read(path):
    with path.open(newline='', encoding='utf-8') as stream:
        return list(csv.DictReader(stream))


def col(rows, name):
    return np.array([float(row[name]) for row in rows])


def vec(rows, prefix, axes='xyz'):
    return np.column_stack([col(rows, prefix+'_'+a) for a in axes])


def rms(values):
    return float(np.sqrt(np.mean(np.asarray(values)**2)))


def finish(fig, name, output):
    for axis in fig.axes:
        if axis.axison:
            axis.grid(alpha=.18)
    fig.savefig(output / name, metadata={'CreationDate':None, 'ModDate':None})
    plt.close(fig)


def rate_errors(nav, measurements):
    imu = {round(float(r['time_s']),6):r for r in measurements if r['sensor']=='imu'}
    values = []
    for row in nav:
        key = round(float(row['time_s']),6)
        if key not in imu:
            continue
        sample = imu[key]
        measured = np.array([float(sample['gyro_body_radps_'+a]) for a in 'xyz'])
        truth = np.array([float(sample['gyro_true_body_radps_'+a]) for a in 'xyz'])
        bias = np.array([float(row['gyro_bias_est_body_radps_'+a]) for a in 'xyz'])
        gain = np.eye(3)+np.diag([float(row['gyro_scale_factor_est_'+a]) for a in 'xyz'])
        for (i,j),name in zip([(0,1),(0,2),(1,0),(1,2),(2,0),(2,1)], ['xy','xz','yx','yz','zx','zy']):
            gain[i,j] = float(row['gyro_misalignment_est_rad_'+name])
        values.append([float(row['time_s']), np.rad2deg(np.linalg.norm(measured-truth)),
            np.rad2deg(np.linalg.norm(np.linalg.solve(gain, measured-bias)-truth))])
    return np.asarray(values)


def model_figures(output):
    fig, ax = plt.subplots(figsize=(7.5,3.1))
    ax.set(xlim=(0,1),ylim=(0,1)); ax.axis('off')
    def box(x,y,w,h,title,body,color):
        ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle='round,pad=0.012',
            facecolor=color+'12',edgecolor=color,linewidth=1.1))
        ax.text(x,y+h*.26,title,ha='center',va='center',fontsize=8.5,fontweight='bold',color=color)
        ax.text(x,y-h*.12,body,ha='center',va='center',fontsize=8.1,linespacing=1.35,color='#203544')
    def arrow(a,b,color=BLUE,style='-',points=()):
        if points:
            v=[a,*points]; ax.plot([p[0] for p in v],[p[1] for p in v],color=color,lw=1,ls=style); a=v[-1]
        ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=9,color=color,lw=1,linestyle=style))
    box(.50,.925,.94,.105,'MISSION AND COMMANDS','mission.py / Events.py / ascent.py / recovery.py',BLUE)
    box(.16,.655,.28,.235,'PHYSICAL MODELS','environment / propulsion / aero\nmass / CG / inertia',GREEN)
    box(.50,.655,.28,.235,'DYNAMICS AND TRUTH','forces / moments / RK4\nposition / velocity / q / rate',BLUE)
    box(.84,.655,.27,.235,'SENSOR GENERATION','IMU / GNSS / LiDAR / star tracker\nerrors / latency / availability',GREEN)
    box(.16,.225,.28,.245,'SEPARATION','conservative state map\nbooster and upper child',BLUE)
    box(.50,.225,.28,.245,'NAVIGATION (MEKF)','gyro quaternion propagation\nGNSS / LiDAR / gated star fixes\n30-state error covariance',ORANGE)
    box(.84,.225,.27,.245,'EVALUATION + OUTPUT','truth / estimates / gyro baseline\nevents / full covariance / errors',BLUE)
    arrow((.16,.86),(.16,.79)); arrow((.50,.86),(.50,.79))
    arrow((.312,.655),(.346,.655)); arrow((.654,.655),(.69,.655))
    ax.add_patch(FancyArrowPatch((.25,.36),(.42,.526),arrowstyle='<->',mutation_scale=9,color=BLUE,lw=1))
    arrow((.80,.526),(.58,.36),GREEN); arrow((.654,.225),(.69,.225),ORANGE)
    arrow((.64,.54),(.84,.089),BLUE,'--',points=((.671,.51),(.671,.03),(.84,.03)))
    ax.text(.75,.035,'truth for comparison',ha='center',va='bottom',fontsize=6.8,color=BLUE)
    fig.tight_layout(pad=.25); finish(fig,'simulation_architecture.pdf',output)
    vehicle=build_falcon9(); lower,_=vehicle.separate(); model=lower.mass_model()
    minimum=model.propellant_bounds(None)[0][0]
    consumed=np.linspace(0,410900-minimum,150)
    props=[model.mass_properties(None,[410900-x,0]) for x in consumed]
    fig,axes=plt.subplots(1,2,figsize=(7.3,2.1),layout='constrained')
    for ax,values,label in [(axes[0],[p.cg_station_body_m[0] for p in props],'CG station from base [m]'),
        (axes[1],[p.inertia_body_kgm2[1,1]/1e6 for p in props],'Transverse inertia [$10^6$ kg m$^2$]')]:
        ax.plot(consumed/1000,values,color=BLUE); ax.set(xlabel='Consumed stage-1 propellant [t]',ylabel=label)
    finish(fig,'mass_properties_depletion.pdf',output)
    lines=(ROOT/'mass_properties.py').read_text().splitlines()
    first=next(i for i,s in enumerate(lines) if '    inertia = np.zeros((3, 3))' in s)
    last=next(i for i in range(first+1,len(lines)) if '    return MassProperties(' in lines[i])
    display=[]
    for k,line in enumerate(lines[first:last]):
        indent=' '*(len(line)-len(line.lstrip())+4)
        parts=textwrap.wrap(line,width=84,subsequent_indent=indent,break_long_words=False,
            break_on_hyphens=False,replace_whitespace=False,drop_whitespace=True) or ['']
        display.extend((str(first+k+1) if j==0 else '>',p) for j,p in enumerate(parts))
    fig=plt.figure(figsize=(7.3,.48+.135*len(display))); fig.patch.set_facecolor('#F3F6F9')
    fig.text(.025,.97,f'mass_properties.py | lines {first+1}-{last}',va='top',weight='bold',color=BLUE,fontsize=10)
    for k,(number,line) in enumerate(display):
        y=1-(.38+k*.135)/fig.get_figheight()
        fig.text(.055,y,number,ha='right',va='top',fontfamily='DejaVu Sans Mono',fontsize=8,color=GRAY)
        fig.text(.069,y,line,va='top',fontfamily='DejaVu Sans Mono',fontsize=8,color='#203544')
    finish(fig,'code_mass_assembly.pdf',output)


def result_figures(run, output):
    truth=read(run/'truth.csv'); upper=read(run/'stage2_truth.csv')
    nav=read(run/'navigation.csv'); ref=read(run/'validation.csv'); measurements=read(run/'measurements.csv')
    time,ut,nt=col(truth,'time_s'),col(upper,'time_s'),col(nav,'time_s')
    p=vec(truth,'cg_position_ecef_m')
    up=p/np.linalg.norm(p,axis=1)[:,None]
    tilt=np.rad2deg(np.arccos(np.clip(np.sum(np.array([body_to_ecef_matrix(q)[:,0] for q in
        vec(truth,'quaternion_body_to_ecef','wxyz')])*up,axis=1),-1,1)))
    fig,ax=plt.subplots(3,2,figsize=(7.3,4.1),layout='constrained')
    ax[0,0].plot(time,col(truth,'cg_altitude_m')/1000,color=GREEN,label='Stack / booster')
    ax[0,0].plot(ut,col(upper,'cg_altitude_m')/1000,color=BLUE,label='Upper + Dragon')
    ax[0,0].set_ylabel('CG altitude [km]'); ax[0,0].legend(fontsize=6.5)
    for rows,t,color in [(truth,time,GREEN),(upper,ut,BLUE)]:
        ax[0,1].plot(t,col(rows,'speed_mps')/1000,color=color)
        ax[1,1].plot(t,col(rows,'thrust_n')/1e6,color=color)
    ax[0,1].set_ylabel('Ground-relative speed [km/s]'); ax[1,1].set_ylabel('Thrust [MN]')
    ax[1,0].plot(time,col(truth,'mass_kg')/1000,color=GREEN,label='Vehicle')
    ax[1,0].plot(time,col(truth,'stage1_propellant_kg')/1000,color=ORANGE,label='Stage-1 propellant')
    ax[1,0].set_ylabel('Mass [t]'); ax[1,0].legend(fontsize=6.5)
    ax[2,0].plot(time,col(truth,'dynamic_pressure_pa')/1000,color=BLUE); ax[2,0].set_ylabel('Dynamic pressure [kPa]')
    ax[2,1].plot(time,np.rad2deg(col(truth,'angle_of_attack_rad')),color=ORANGE,label='Flow angle')
    ax[2,1].plot(time,tilt,color=GREEN,label='Nose tilt from up'); ax[2,1].set_ylabel('Angle [deg]'); ax[2,1].legend(fontsize=6.5)
    for a in ax.flat: a.set_xlim(0,600)
    for a in ax[-1]: a.set_xlabel('Time [s]')
    finish(fig,'results_flight_profiles.pdf',output)
    fig,ax=plt.subplots(1,2,figsize=(7.3,2.0),layout='constrained')
    for rows,color,label in [(truth,GREEN,'Stack / booster'),(upper,BLUE,'Upper + Dragon')]:
        position=vec(rows,'cg_position_ecef_m'); radius=np.linalg.norm(position,axis=1)
        downrange=6371000*np.arccos(np.clip(position[:,2]/radius,-1,1))
        ax[0].plot(-position[:,1]/1000,position[:,0]/1000,color=color,label=label)
        ax[1].plot(downrange/1000,(radius-6371000)/1000,color=color,label=label)
    ax[0].scatter([0,6],[0,-14],s=12,color=[GRAY,ORANGE]); ax[0].set(xlabel='Launch-frame east [km]',ylabel='Launch-frame north [km]')
    ax[1].set(xlabel='Great-circle downrange [km]',ylabel='CG altitude [km]'); ax[1].legend(fontsize=7)
    finish(fig,'results_trajectories.pdf',output)
    fig,ax=plt.subplots(3,1,figsize=(7.3,4.2),layout='constrained',sharex=True)
    for rows,color,label in [(truth,GREEN,'Stack / booster'),(upper,BLUE,'Upper + Dragon')]:
        t=col(rows,'time_s'); position=vec(rows,'cg_position_ecef_m')
        unit=position/np.linalg.norm(position,axis=1)[:,None]
        ax[0].plot(t,6371*np.arccos(np.clip(unit[:,2],-1,1)),color=color,label=label)
        if 'cg_velocity_ecef_mps_x' in rows[0]:
            velocity=vec(rows,'cg_velocity_ecef_mps'); radial=np.sum(velocity*unit,axis=1)
            tangent=np.linalg.norm(velocity-radial[:,None]*unit,axis=1)
            gamma=np.rad2deg(np.arctan2(radial,tangent)); gamma[np.linalg.norm(velocity,axis=1)<1]=np.nan
            ax[1].plot(t,gamma,color=color)
    omega=np.rad2deg(vec(truth,'omega_body_radps'))
    for j,color in enumerate([BLUE,ORANGE,GREEN]): ax[2].plot(time,omega[:,j],color=color,lw=.9,label='$\\omega_'+ 'xyz'[j]+'$')
    ax[2].plot(time,np.linalg.norm(omega,axis=1),color=GRAY,lw=.8,ls='--',label='Magnitude')
    ax[0].set_ylabel('Downrange [km]'); ax[0].legend(fontsize=7,ncol=2)
    ax[1].set_ylabel('Flight-path angle [deg]'); ax[2].set(ylabel='Body rate [deg/s]',xlabel='Time [s]',xlim=(0,600))
    ax[2].legend(fontsize=7,ncol=4)
    for a in ax: a.axvline(144,color=GRAY,ls=':',lw=.7)
    finish(fig,'results_kinematics.pdf',output)
    fig,ax=plt.subplots(1,2,figsize=(7.3,2.15),layout='constrained')
    for a,field,scale,label in [(ax[0],'altitude_m',1e-3,'CG altitude [km]'),(ax[1],'speed_mps',1,'Speed [m/s]')]:
        a.plot(col(ref,'time_s'),col(ref,'simulated_'+field)*scale,color=BLUE,label='Simulation')
        a.plot(col(ref,'time_s'),col(ref,'reference_'+field)*scale,'--',color=ORANGE,label='Video readout (E)')
        a.set(xlabel='Time [s]',ylabel=label,xlim=(0,140)); a.legend(fontsize=7)
    finish(fig,'results_ascent_reference.pdf',output)
    pos=vec(nav,'position_error_ecef_m'); vel=vec(nav,'velocity_error_ecef_mps')
    sigma=np.column_stack([col(nav,f'sigma_{j}') for j in range(3)])
    rates=rate_errors(nav,measurements); star=read(run/'startracker_innovations.csv')
    fig,ax=plt.subplots(3,2,figsize=(7.3,4.7),layout='constrained')
    ax[0,0].plot(nt,np.linalg.norm(pos,axis=1),color=BLUE); ax[0,0].set_ylabel('3D position error [m]')
    ax[0,1].plot(nt,np.linalg.norm(vel,axis=1),color=BLUE); ax[0,1].set_ylabel('3D velocity error [m/s]')
    ax[1,0].plot(nt,col(nav,'gyro_only_attitude_error_deg'),color=GRAY,label='Gyro only')
    ax[1,0].plot(nt,col(nav,'attitude_error_deg'),color=BLUE,label='MEKF'); ax[1,0].set_ylabel('Attitude error [deg]'); ax[1,0].legend(fontsize=7)
    ax[1,1].plot(rates[:,0],rates[:,1],color=GRAY,lw=.6,label='Raw')
    ax[1,1].plot(rates[:,0],rates[:,2],color=GREEN,lw=.5,label='Calibrated'); ax[1,1].set_ylabel('Rate error [deg/s]'); ax[1,1].legend(fontsize=7)
    for j,color in enumerate([BLUE,ORANGE,GREEN]): ax[2,0].plot(nt,pos[:,j]/sigma[:,j],color=color,lw=.7,label='xyz'[j])
    ax[2,0].axhspan(-3,3,color=GRAY,alpha=.15); ax[2,0].set_ylabel('Position error / sigma'); ax[2,0].legend(fontsize=7,ncol=3)
    ax[2,1].plot(nt,col(nav,'gnss_nis'),color=GRAY,lw=.55,label='GNSS (6D)')
    ax[2,1].scatter(col(star,'time_s'),col(star,'nis'),s=3,color=BLUE,label='Star (3D)')
    ax[2,1].axhline(16.27,color=ORANGE,lw=.7,ls='--',label='Star gate'); ax[2,1].set_ylabel('NIS'); ax[2,1].set_yscale('symlog',linthresh=1)
    ax[2,1].legend(fontsize=6,ncol=2)
    ax[2,1].set_ylim(bottom=0)
    for a in ax.flat: a.set_xlim(0,nt[-1])
    for a in ax[-1]: a.set_xlabel('Time [s]')
    finish(fig,'results_navigation.pdf',output)
    fig,ax=plt.subplots(3,1,figsize=(7.3,4.35),layout='constrained',sharex=True)
    ax[0].plot(time,col(truth,'cg_station_m'),color=GREEN); ax[0].set_ylabel('CG station [m]')
    for name,label,color in [('inertia_kgm2_x','$I_{xx}$',GREEN),('inertia_kgm2_y','$I_{yy}=I_{zz}$',BLUE)]:
        ax[1].plot(time,col(truth,name),color=color,label=label)
    ax[1].set(yscale='log',ylabel='CG inertia [kg m$^2$]'); ax[1].legend(fontsize=7)
    counts={}
    for y,sensor in enumerate(['gnss','lidar','startracker']):
        rows=[r for r in measurements if r['sensor']==sensor]; t=col(rows,'time_s')
        valid=np.array([r['valid'].lower()=='true' for r in rows]); counts[sensor]={'emitted':len(rows),'valid':int(valid.sum())}
        ax[2].scatter(t[~valid],np.full((~valid).sum(),y),s=1.2,color='#d0d4da')
        ax[2].scatter(t[valid],np.full(valid.sum(),y),s=2,color=GREEN)
    ax[2].set(yticks=[0,1,2],yticklabels=['GNSS','LiDAR','Star tracker'],ylim=(-.5,2.7),xlabel='Time [s]')
    accepted=int(col(nav,'startracker_accepted')[-1]); rejected=int(col(nav,'startracker_rejected')[-1])
    ax[2].text(.99,.93,f'Star tracker: {accepted} accepted; {rejected} rejected',transform=ax[2].transAxes,ha='right',va='top',fontsize=8)
    for a in ax: a.axvline(144,color=GRAY,ls='--',lw=.7); a.set_xlim(0,time[-1])
    finish(fig,'appendix_properties_availability.pdf',output)
    return truth,upper,nav,measurements,counts,rates


def validation_figure(output):
    data=json.loads((REPORT/'navigation_validation.json').read_text())
    g=data['gyro_calibration']; trace=np.array(g['trace']); c=data['covariance_consistency']
    fig,ax=plt.subplots(1,3,figsize=(7.3,2.05),layout='constrained')
    ax[0].plot(trace[:,0],trace[:,1],color=BLUE); ax[0].set(xlabel='Time [s]',ylabel='Attitude error [deg]')
    ax[1].semilogy(trace[:,0],trace[:,3],color=BLUE,label='Gyro scale')
    ax[1].semilogy(trace[:,0],trace[:,4],color=GREEN,label='Cross-axis [rad]')
    ax[1].set(xlabel='Time [s]',ylabel='Calibration error norm'); ax[1].legend(fontsize=6)
    ax[2].scatter([r['seed'] for r in c['per_seed']],[r['position_nees'] for r in c['per_seed']],s=10,color=BLUE)
    ax[2].axhline(c['position_anees'],color=GREEN,label='Mean')
    ax[2].set(xlabel='Independent trial',ylabel='Final position NEES'); ax[2].legend(fontsize=6)
    finish(fig,'results_navigation_validation.pdf',output)


def metrics(run, baseline, coarse, loaded):
    truth,upper,nav,measurements,counts,rates=loaded
    summary=json.loads((run/'summary.json').read_text()); events=json.loads((run/'events.json').read_text())
    cov=np.load(run/'navigation_covariance.npz'); matrices=cov['covariance']
    assert matrices.shape==(len(nav),30,30) and np.isfinite(matrices).all()
    assert np.allclose(cov['time_s'],col(nav,'time_s'),atol=1e-8,rtol=0)
    symmetry=float(np.max(np.abs(matrices-matrices.transpose(0,2,1))))
    eigen=float(np.min(np.linalg.eigvalsh(matrices)))
    assert symmetry<1e-10 and eigen>=-1e-12
    pos=vec(nav,'position_error_ecef_m'); vel=vec(nav,'velocity_error_ecef_mps')
    att=col(nav,'attitude_error_deg'); sigma=np.column_stack([col(nav,f'sigma_{j}') for j in range(3)])
    star=read(run/'startracker_innovations.csv')
    data={'run':str(run.relative_to(ROOT)).replace('\\','/'),'runtime_revision':summary['git_commit'],
        'runtime_hashes_verified':len(summary['input_sha256']), 'position_rms_m':rms(np.linalg.norm(pos,axis=1)),
        'velocity_rms_mps':rms(np.linalg.norm(vel,axis=1)), 'attitude_rms_deg':rms(att),
        'attitude_peak_deg':float(att.max()),'attitude_peak_s':float(col(nav,'time_s')[att.argmax()]),
        'attitude_final_deg':float(att[-1]),'gyro_only_final_deg':float(col(nav,'gyro_only_attitude_error_deg')[-1]),
        'raw_rate_rms_degps':rms(rates[:,1]),'calibrated_rate_rms_degps':rms(rates[:,2]),'rate_sample_count':len(rates),
        'position_3sigma_time_fraction':np.mean(np.abs(pos)<=3*sigma,axis=0).tolist(),
        'star_accepted':int(col(nav,'startracker_accepted')[-1]),'star_rejected':int(col(nav,'startracker_rejected')[-1]),
        'star_nis_mean':float(np.mean(col(star,'nis'))),'star_nis_max':float(np.max(col(star,'nis'))),
        'gnss_nis_mean':float(np.nanmean(col(nav,'gnss_nis'))),'covariance_shape':list(matrices.shape),
        'covariance_min_eigenvalue':eigen,'covariance_max_asymmetry':symmetry,
        'landing':next(e['landing_metrics'] for e in events if 'landing_metrics' in e),
        'reference_comparison':summary['reference_comparison'],'burns':summary['recovery_burns'],
        'sensor_counts':counts,'events':[{'event':e['event'],'time_s':e['time_s']} for e in events],
        'input_sha256':summary['input_sha256'], 'source_hashes':{f:digest(run/f) for f in FILES if (run/f).exists()}}
    if baseline and baseline.exists():
        old=read(baseline/'navigation.csv')
        data['paired_baseline']={'run':str(baseline.relative_to(ROOT)).replace('\\','/'),
            'position_rms_m':rms(np.linalg.norm(vec(old,'position_error_ecef_m'),axis=1)),
            'velocity_rms_mps':rms(np.linalg.norm(vec(old,'velocity_error_ecef_mps'),axis=1)),
            'attitude_rms_deg':rms(col(old,'attitude_error_deg')),
            'truth_byte_identical':digest(run/'truth.csv')==digest(baseline/'truth.csv'),
            'upper_truth_byte_identical':digest(run/'stage2_truth.csv')==digest(baseline/'stage2_truth.csv')}
        previous=read(baseline/'measurements.csv')
        keys=[k for k in previous[0] if k not in {'measurement_epoch_s','latency_s'}]
        data['paired_baseline']['common_measurements_identical']=len(previous)==len(measurements) and all(
            all(a[k]==b[k] for k in keys) for a,b in zip(previous,measurements))
        assert all(data['paired_baseline'][k] for k in ['truth_byte_identical','upper_truth_byte_identical','common_measurements_identical'])
    if coarse and coarse.exists():
        data['coarse']={'run':str(coarse.relative_to(ROOT)).replace('\\','/'),
            'landing':next(e['landing_metrics'] for e in json.loads((coarse/'events.json').read_text()) if 'landing_metrics' in e)}
    return data


def export_numbers(data):
    values={'PositionRMS':('position_rms_m','.3f'),'VelocityRMS':('velocity_rms_mps','.4f'),
        'AttitudeRMS':('attitude_rms_deg','.3f'),'AttitudeFinal':('attitude_final_deg','.5f'),
        'AttitudePeak':('attitude_peak_deg','.3f'),'AttitudePeakTime':('attitude_peak_s','.1f'),
        'GyroOnlyFinal':('gyro_only_final_deg','.3f'),'RawRateRMS':('raw_rate_rms_degps','.5f'),
        'CalibratedRateRMS':('calibrated_rate_rms_degps','.5f'),'GnssNIS':('gnss_nis_mean','.2f'),
        'StarNIS':('star_nis_mean','.2f'),'StarMaxNIS':('star_nis_max','.2f')}
    lines=['\\newcommand{\\'+name+'}{'+format(data[key],fmt)+'}' for name,(key,fmt) in values.items()]
    for name,key in [('StarAccepted','star_accepted'),('StarRejected','star_rejected'),('RateSamples','rate_sample_count')]:
        lines.append('\\newcommand{\\'+name+'}{'+str(data[key])+'}')
    for axis,value in zip('XYZ',data['position_3sigma_time_fraction']):
        lines.append(r'\newcommand{\Coverage'+axis+'}{'+f'{value*100:.1f}'+r'\%}')
    (REPORT/'content/results_numbers.tex').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')


def pack(paths, archive):
    records={}
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as output:
        for label,path in paths.items():
            for name in FILES:
                source=path/name
                if not source.exists(): continue
                key=label+'/'+name; output.write(source,key); records[key]=digest(source)
    return {'file':archive.name,'sha256':digest(archive),'records':records}


def unpack(archive, destination, ledger):
    if digest(archive)!=ledger['sha256']: raise ValueError('Result archive checksum mismatch.')
    destination=destination.resolve(); destination.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(archive) as stream:
        if set(stream.namelist())!=set(ledger['records']): raise ValueError('Unexpected archive contents.')
        for member in stream.infolist():
            target=(destination/member.filename).resolve()
            if not target.is_relative_to(destination): raise ValueError('Unsafe archive path.')
            content=stream.read(member)
            if hashlib.sha256(content).hexdigest()!=ledger['records'][member.filename]: raise ValueError('Record checksum mismatch.')
            if target.exists() and target.read_bytes()!=content: raise FileExistsError(target)
            target.parent.mkdir(parents=True,exist_ok=True); target.write_bytes(content)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run',type=Path,default=ROOT/'results/final-submission-20261110-nav30')
    parser.add_argument('--baseline',type=Path,default=ROOT/'results/gyro-final-review-20261110')
    parser.add_argument('--coarse',type=Path,default=ROOT/'results/two-stage-final-check')
    parser.add_argument('--output-dir',type=Path,default=REPORT/'images')
    parser.add_argument('--pack',action='store_true'); parser.add_argument('--unpack',action='store_true')
    parser.add_argument('--validate',action='store_true')
    args=parser.parse_args(); manifest=json.loads((REPORT/'figure_sources.json').read_text())
    archive=REPORT/'reproduction_data.zip'
    if args.unpack:
        base=ROOT/'results/submission-reproduction'
        unpack(archive,base,manifest['reproduction_archive'])
        args.run, args.baseline,args.coarse=[base/k for k in ['nominal','baseline','coarse']]
    if args.validate:
        from validate_navigation import run_all
        validation=run_all()
        if not all(validation['checks'].values()): raise ValueError('Independent validation failed.')
        (REPORT/'navigation_validation.json').write_text(json.dumps(validation,indent=2)+'\n',encoding='utf-8',newline='\n')
    canonical,changed_endings=verify_inputs(args.run,manifest.get('runtime_text_lf_sha256',{}))
    args.output_dir.mkdir(parents=True,exist_ok=True)
    loaded=result_figures(args.run,args.output_dir)
    data=metrics(args.run,args.baseline,args.coarse,loaded)
    if args.unpack:
        data['run']=manifest['integrated_results']['run']
        for key in ['paired_baseline','coarse']:
            data[key]['run']=manifest['integrated_results'][key]['run']
    model_figures(args.output_dir); validation_figure(args.output_dir); export_numbers(data)
    manifest['source_revision']=data['runtime_revision']; manifest['integrated_results']=data
    manifest['runtime_text_lf_sha256']=canonical
    manifest.pop('navigation_review',None)
    manifest['model_input_sha256']={name:digest(ROOT/name) for name in
        ['vehicle_parameters.py','vehicle_mass_model.py','mass_properties.py','data/vehicle_parameters.json']}
    manifest['model_figures'][1]['basis']='Implemented interfaces and nominal 30-component MEKF error state.'
    manifest['result_figures']=['results_flight_profiles.pdf','results_trajectories.pdf','results_kinematics.pdf',
        'results_ascent_reference.pdf','results_navigation.pdf','results_navigation_validation.pdf']
    manifest['final_appendix']={'nominal_run':data['run'],'runtime_hashes_verified':data['runtime_hashes_verified'],
        'sensor_counts':data['sensor_counts'],'source_files':data['source_hashes'],
        'code_excerpts':[{'file':name,'first_line':first,'last_line':last,'sha256':digest(ROOT/name)}
            for name,first,last in [('load_interface.py',27,41),('dynamics.py',163,176),('VehicleData.py',168,184)]]}
    manifest['navigation_validation_sha256']=digest(REPORT/'navigation_validation.json')
    manifest['reproducer_sha256']=digest(Path(__file__))
    if args.pack: manifest['reproduction_archive']=pack({'nominal':args.run,'baseline':args.baseline,'coarse':args.coarse},archive)
    manifest['generated_figure_sha256']={p.name:digest(p) for p in args.output_dir.glob('*.pdf') if p.name in
        manifest['result_figures']+['simulation_architecture.pdf','mass_properties_depletion.pdf','code_mass_assembly.pdf','appendix_properties_availability.pdf']}
    (REPORT/'figure_sources.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:v for k,v in data.items() if k not in {'input_sha256','source_hashes','events'}},indent=2))


if __name__=='__main__':
    main()
