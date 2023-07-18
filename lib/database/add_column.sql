ALTER TABLE simulation_param
  ADD "division_time_hours" numeric;

ALTER TABLE simulation_param
  ADD "p_division" numeric;


ALTER table simulation_param
ADD COLUMN id UUID not null DEFAULT uuid_generate_v4();

-- Install the uuid-ossp extension (if not already installed)
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";

-- Enable the extension
SELECT pg_catalog.pg_extension.* FROM pg_catalog.pg_extension WHERE extname = 'uuid-ossp';
