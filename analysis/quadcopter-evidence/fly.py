# One flight on real PX4 in SIH: climb, then a roll step and a roll-rate step.
import glob, json, math, os, shutil, subprocess, sys, time
from pymavlink import mavutil
tag = sys.argv[1]; cfg = json.load(open("/out/plant.json")); root = "/src/build/px4_sitl_nodds"; rootfs = root + "/rootfs"
for f in glob.glob(rootfs + "/parameters*.bson") + glob.glob(rootfs + "/dataman"): os.remove(f)
shutil.rmtree(rootfs + "/log", ignore_errors=True)
params = dict(SIH_MASS=cfg["mass"], SIH_IXX=cfg["ixx"], SIH_IYY=cfg["iyy"], SIH_IZZ=cfg["izz"], SIH_T_MAX=cfg["tmax"], SIH_Q_MAX=cfg["qmax"], SIH_L_ROLL=cfg["lroll"], SIH_L_PITCH=cfg["lroll"], SIH_T_TAU=cfg["tau"], SIH_KDV=0.3, SIH_KDW=0,
    SDLOG_PROFILE=25, COM_RCL_EXCEPT=4, NAV_DLL_ACT=0, COM_DISARM_PRFLT=-1, MPC_THR_HOVER=round(min(.8, max(.1, cfg["hover"])), 3), COM_OF_LOSS_T=5, **cfg["gains"][tag])
env = dict(os.environ, PX4_SYS_AUTOSTART="10040", PX4_SIMULATOR="sihsim", PX4_SIM_MODEL="quadx", HEADLESS="1", **{"PX4_PARAM_" + k: str(v) for k, v in params.items()})
out = open(f"/out/{tag}.px4.log", "w"); px4 = subprocess.Popen([root + "/bin/px4", "-d"], env=env, stdout=out, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
m = mavutil.mavlink_connection("udpin:127.0.0.1:14540", source_system=255); m.wait_heartbeat(timeout=60); print("heartbeat", flush=True)
MODE, ARM = mavutil.mavlink.MAV_CMD_DO_SET_MODE, mavutil.mavlink.MAV_CMD_COMPONENT_ARM_DISARM
alt = 30.0
def hold(): m.mav.set_position_target_local_ned_send(0, 1, 1, mavutil.mavlink.MAV_FRAME_LOCAL_NED, 0b110111111000, 0, 0, -alt, 0, 0, 0, 0, 0, 0, 0, 0)
def att(r): m.mav.set_attitude_target_send(0, 1, 1, 0b00000111, [math.cos(r / 2), math.sin(r / 2), 0, 0], 0, 0, 0, cfg["hover"])
def rate(p): m.mav.set_attitude_target_send(0, 1, 1, 0b10000000, [1, 0, 0, 0], p, 0, 0, cfg["hover"])
def stream(fn, seconds):
    end = time.time() + seconds
    while time.time() < end:
        fn(); time.sleep(.01)
        while m.recv_match(blocking=False): pass
armed = lambda: bool(m.messages.get("HEARTBEAT") and m.messages["HEARTBEAT"].base_mode & mavutil.mavlink.MAV_MODE_FLAG_SAFETY_ARMED)
start = time.time()
while not armed():
    if time.time() - start > 90: print("never armed"); px4.terminate(); sys.exit(2)
    stream(hold, 1)
    m.mav.command_long_send(1, 1, MODE, 0, 1, 6, 0, 0, 0, 0, 0); m.mav.command_long_send(1, 1, ARM, 0, 1, 0, 0, 0, 0, 0, 0)
print("armed after", round(time.time() - start, 1), "s", flush=True)
stream(hold, 24); pos = m.messages.get("LOCAL_POSITION_NED"); print("at altitude", pos and round(-pos.z, 2), flush=True)
step = math.radians(cfg["step"])
stream(lambda: att(0), 2); stream(lambda: att(step), 3); stream(lambda: att(0), 3); stream(hold, 6)
# 24 rate steps of each sign in a pattern that keeps the aircraft near level, to average the noise away.
stream(lambda: rate(0), 1)
for sign in [1, -1, -1, 1] * 6:
    stream(lambda: rate(sign * math.radians(30)), .4); stream(lambda: rate(0), .8)
stream(hold, 6)
# A logarithmic sweep of roll rate at about 1 degree of attitude, 0.7 to 16 Hz, to identify the loop itself.
sweep, f0, f1 = 28.0, .7, 16.0; k = math.log(f1 / f0) / sweep; begin = time.time()
def chirp():
    x = time.time() - begin; rate(math.radians(10) * max(1, f0 * math.exp(k * x) / 1.5) * math.sin(2 * math.pi * f0 * (math.exp(k * x) - 1) / k))
stream(lambda: rate(0), .5); begin = time.time(); stream(chirp, sweep); stream(lambda: rate(0), .5); stream(hold, 4)
m.mav.command_long_send(1, 1, ARM, 0, 0, 21196, 0, 0, 0, 0, 0); stream(hold, 2); px4.terminate(); px4.wait(timeout=20)
logs = sorted(glob.glob(rootfs + "/log/**/*.ulg", recursive=True)); print("logs", logs, flush=True); shutil.copy(logs[-1], f"/out/{tag}.ulg")
from pyulog import ULog
u = ULog(f"/out/{tag}.ulg"); print("topics", len(u.data_list), flush=True)
def series(name, fields, multi=0):
    for d in u.data_list:
        if d.name == name and d.multi_id == multi:
            t = d.data["timestamp"] / 1e6
            return {"t": t.tolist(), **{k: d.data[k].tolist() for k in fields if k in d.data}}
    return None
Q = [f"q[{i}]" for i in range(4)]
res = dict(att=series("vehicle_attitude", Q), att_sp=series("vehicle_attitude_setpoint", ["q_d[0]", "q_d[1]", "q_d[2]", "q_d[3]", "thrust_body[2]"]), rate=series("vehicle_angular_velocity", ["xyz[0]", "xyz[1]", "xyz[2]"]), rate_sp=series("vehicle_rates_setpoint", ["roll", "pitch", "yaw", "thrust_body[2]"]), motors=series("actuator_motors", [f"control[{i}]" for i in range(4)]), torque=series("vehicle_torque_setpoint", ["xyz[0]", "xyz[1]", "xyz[2]"]), pos=series("vehicle_local_position", ["z", "vz"]))
json.dump(res, open(f"/out/{tag}.json", "w")); print({k: (len(v["t"]) if v else None) for k, v in res.items()}, flush=True)
