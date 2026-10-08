// Tests setup
pub static TEST_TIME: f64 = 3.15576e17;
// 9 Gyr: inside the 2026 STAREVOL solar track, which ends at 9.67 Gyr (before TEST_TIME).
pub static TEST_TIME_mode_fromfile: f64 = 2.840184e17;
pub static TEST_DISK_LIFETIME: f64 = 78325963200000.;
pub static DISK_IS_DISSIPATED: bool = !(TEST_TIME < TEST_DISK_LIFETIME);
