-- This script is executed when the Docker container is first created.

-- Create the table for driver Thiele/Small parameters
CREATE TABLE drivers (
    driver_id VARCHAR(255) PRIMARY KEY,
    manufacturer VARCHAR(255),
    model_name VARCHAR(255),
    fs_hz REAL,
    qts REAL,
    vas_liters REAL,
    sd_sq_meters REAL,
    re_ohms REAL,
    le_mh REAL,
    xmax_mm REAL
);

-- Insert a sample driver for testing purposes
INSERT INTO drivers (
    driver_id, manufacturer, model_name, fs_hz, qts, vas_liters, sd_sq_meters, re_ohms, le_mh, xmax_mm
) VALUES (
    'test_driver_001',
    'Test Audio',
    'Model-X',
    35.5,
    0.35,
    140.0,
    0.053,
    5.8,
    0.4,
    8.5
); 