#!/usr/bin/python3
print("================================")
print("      DATA RATE PER UNIT        ")

# ASSUME FLOAT32 FOR ALL MEASUREMENTS
float_size = 32; # float 32 size
int_size = 16;
bool_size = 8;
protocol_overhead = 2; # MIL-STD-1553 protocol overhead
# A factor of 2 was conservatively adopted to account for 
# MIL-STD-1553 framing, synchronization, command/status words, 
# and bus-management overhead.

# ==== STAR TRACKER ==== #
st_update_rate = 10; # [Hz]
st_attitude_packet = 4 * float_size; # assume quaternions
st_angular_rates_packet = 3 * float_size; # 3-axis
st_housekeeping_packet = float_size; # advanced component so assume secondary outputs

st_bitrate = (st_attitude_packet + st_angular_rates_packet + st_housekeeping_packet) * protocol_overhead * st_update_rate;
st_n = 2;
st_bitrate_tot = st_bitrate * st_n;
print(f"STAR TRACKER    : {int(st_bitrate)} bps");

# ==== GYRO ==== #
g_axes = 2; # how many axes per gyroscope
g_update_rate = 75;
#g_update_rate_lo = 50; # TODO update rate [Hz]
#g_update_rate_hi = 100; # TODO justify BOHHH

g_bitrate = g_update_rate * g_axes * float_size * protocol_overhead;
g_n = 3;
g_bitrate_tot = g_bitrate * g_n;
print(f"GYRO            : {int(g_bitrate)} bps");

#g_bitrate_lo = g_update_rate_lo * g_axes * float_size * protocol_overhead;
#g_bitrate_hi = g_update_rate_hi * g_axes * float_size * protocol_overhead;
#print(f"GYRO (LO)       : {g_bitrate_lo} bps");
#print(f"GYRO (HI)       : {g_bitrate_hi} bps");

# ==== COARSE SUN SENSOR ==== #
css_channels = 4; # 4 channels reducing to 2 axes
css_update_rate = 2; # TODO BOH

css_bitrate = (int_size) * css_channels * css_update_rate;
css_n = 8;
css_bitrate_tot = css_bitrate * css_n;
print(f"CSS             : {int(css_bitrate)} bps");

# ==== MAG ==== #
mag_update_rate = 2;
mag_axes = 3;

mag_bitrate = mag_axes * mag_update_rate * float_size * protocol_overhead;
mag_n = 2;
mag_bitrate_tot = mag_bitrate * mag_n;
print(f"MAG             : {int(mag_bitrate)} bps");

# ==== GPS ==== # TODO tutto inferred dalla cit. di topstar3000
gps_axes = 3; # X, Y, Z for position and velocity
gps_position_packet = (float_size*2) * gps_axes;
gps_velocity_packet = float_size * gps_axes;
gps_time_packet = float_size * 3; # clock bias, clock drift, GPS time lag

gps_signal_processing_packet = 50;

gps_raw_sats_tracked = 10;
gps_raw_packet_size_per_sat = 5 * float_size;
gps_raw_data_packet = gps_raw_sats_tracked * gps_raw_packet_size_per_sat; # 6-12 sats tracked, est. 5 data segments per satellite

gps_update_rate = 1; # [Hz]

gps_bitrate = ((gps_position_packet + gps_velocity_packet + gps_time_packet) + \
                (gps_signal_processing_packet) + \
                (gps_raw_data_packet)) * protocol_overhead * gps_update_rate;
gps_n = 1;
gps_bitrate_tot = gps_bitrate * gps_n;
print(f"GPS             : {int(gps_bitrate)} bps");

#gps_bitrate_lo = ((gps_position_packet + gps_velocity_packet + gps_time_packet) + \
#                (gps_signal_processing_packet)) * protocol_overhead * gps_update_rate;
#gps_bitrate_hi = ((gps_position_packet + gps_velocity_packet + gps_time_packet) + \
#                (gps_signal_processing_packet) + \
#                (gps_raw_data_packet)) * protocol_overhead * gps_update_rate;
#print(f"GPS (LO)        : {gps_bitrate_lo} bps");
#print(f"GPS (HI)        : {gps_bitrate_hi} bps");

# ==== RW ==== #
rw_out_satur_flag = bool_size;
rw_out_speed = float_size;
rw_in_command = float_size;
#rw_update_rate_lo = 10;
#rw_update_rate_hi = 50;
rw_update_rate = 25;

rw_bitrate = (rw_out_satur_flag + rw_out_speed + rw_in_command) * protocol_overhead * rw_update_rate
rw_n = 4;
rw_bitrate_tot = rw_bitrate * rw_n;
print(f"RW              : {int(rw_bitrate)} bps");

#rw_bitrate_lo = (rw_out_satur_flag + rw_out_speed + rw_in_command) * protocol_overhead * rw_update_rate_lo
#rw_bitrate_lo = (rw_out_satur_flag + rw_out_speed + rw_in_command) * protocol_overhead * rw_update_rate_hi
#print(f"RW (LO)         : {gps_bitrate_lo} bps");
#print(f"RW (HI)         : {gps_bitrate_hi} bps");

# ==== MTQ ==== #
mtq_in_command = float_size;
mtq_update_rate = 3;

mtq_bitrate = mtq_in_command * mtq_update_rate * protocol_overhead;
mtq_n = 3;
mtq_bitrate_tot = mtq_bitrate * mtq_n;
print(f"MTQ             : {int(mtq_bitrate)} bps");

# ==== THRUSTER ==== #
thr_in_command = bool_size;

thr_bitrate = bool_size * protocol_overhead;
thr_n = 4;
thr_bitrate_tot = thr_bitrate * thr_n;
print(f"THRUSTER        : {int(thr_bitrate)} bps");

print("--------------------------------")
#tot_bitrate_lo = st_bitrate + g_bitrate_lo + css_bitrate + mag_bitrate + gps_bitrate_lo;
#tot_bitrate_hi = st_bitrate + g_bitrate_hi + css_bitrate + mag_bitrate + gps_bitrate_hi;
tot_bitrate = st_bitrate + g_bitrate + css_bitrate + mag_bitrate + gps_bitrate + rw_bitrate + mtq_bitrate + thr_bitrate;
print(f"TOTAL PER ITEM  : {int(tot_bitrate)} bps");
#print(f"TOTAL (LO)      : {tot_bitrate_lo} bps");
#print(f"TOTAL (HI)      : {tot_bitrate_hi} bps");
print("--------------------------------")

print("================================")
print("       DATA RATE OVERALL        ")
print(f"STAR TRACKER x {st_n}   : {int(st_bitrate_tot)} bps");
print(f"GYRO x {g_n}           : {int(g_bitrate_tot)} bps");
print(f"CSS x {css_n}            : {int(css_bitrate_tot)} bps");
print(f"MAG x {mag_n}            : {int(mag_bitrate_tot)} bps");
print(f"GPS x {gps_n}            : {int(gps_bitrate_tot)} bps");
print(f"RW x {rw_n}             : {int(rw_bitrate_tot)} bps");
print(f"MTQ x {mtq_n}            : {int(mtq_bitrate_tot)} bps");
print(f"THRUSTER x {thr_n}       : {int(thr_bitrate_tot)} bps");

print("--------------------------------")
tot_bitrate_tot = st_bitrate_tot + g_bitrate_tot + css_bitrate_tot + mag_bitrate_tot + gps_bitrate_tot + rw_bitrate_tot + mtq_bitrate_tot + thr_bitrate_tot;
print(f"OVERALL TOTAL      : {int(tot_bitrate_tot)} bps");
print("--------------------------------")

# ==== MODES ==== #
shm_bitrate = css_bitrate_tot + mag_bitrate_tot + (rw_bitrate_tot/2) + mtq_bitrate_tot + gps_bitrate_tot;
#print(shm_bitrate)
sam_bitrate = css_bitrate_tot + mag_bitrate_tot + rw_bitrate_tot + mtq_bitrate_tot + gps_bitrate_tot + st_bitrate_tot + g_bitrate_tot;
#print(sam_bitrate)
nam_bitrate = rw_bitrate_tot + mtq_bitrate_tot + gps_bitrate_tot + st_bitrate_tot + g_bitrate_tot;
#print(nam_bitrate)
ocm2_bitrate = rw_bitrate_tot + mtq_bitrate_tot + gps_bitrate_tot + st_bitrate_tot + g_bitrate_tot + (thr_bitrate_tot/2);
#print(ocm2_bitrate)
ocm4_bitrate = rw_bitrate_tot + mtq_bitrate_tot + gps_bitrate_tot + st_bitrate_tot + g_bitrate_tot + thr_bitrate_tot;
#print(ocm4_bitrate)
print("================================")
print("             MODES              ")
print(f"SAFE HOLD          : {int(shm_bitrate)} bps");
print(f"STAR ACQUISITION   : {int(sam_bitrate)} bps");
print(f"NORMAL AUTONOMOUS  : {int(nam_bitrate)} bps");
print(f"ORBIT CORRECTION 2 : {int(ocm2_bitrate)} bps");
print(f"ORBIT CORRECTION 4 : {int(ocm4_bitrate)} bps");


print("================================")
