create table if not exists public.pattern_class_output(
    "model_param_id" text not null,
    "simulation_param_uuid" uuid,
    "ssID" int not null,
    "pattern_class" text not null,

    primary key ("model_param_id", "simulation_param_uuid", "ssID")
);

--unique constraint for primary keys
ALTER TABLE pattern_class_output
ADD CONSTRAINT unique_pattern_class_output UNIQUE ("model_param_id", "simulation_param_uuid", "ssID");
--foreign keys to tables
ALTER TABLE pattern_class_output
ADD CONSTRAINT fk_model_param_id
FOREIGN KEY (model_param_id)
REFERENCES model_param (model_param_id);

ALTER TABLE pattern_class_output
ADD CONSTRAINT fk_simulation_param_uuid
FOREIGN KEY (simulation_param_uuid)
REFERENCES simulation_param (simulation_param_uuid);

ALTER TABLE pattern_class_output
ADD CONSTRAINT fk_pattern_class
FOREIGN KEY (pattern_class)
REFERENCES pattern_class_param (pattern_class);

alter table pattern_class_output rename column pattern_class to pattern_class_nogrowth;
alter table pattern_class_output add pattern_class_openboundary text;
alter table pattern_class_output add pattern_class_edgegrowth2 text;

ALTER TABLE pattern_class_output
ADD CONSTRAINT fk_pattern_class_openboundary
FOREIGN KEY (pattern_class_openboundary)
REFERENCES pattern_class_param (pattern_class);

ALTER TABLE pattern_class_output
ADD CONSTRAINT fk_pattern_class_edgegrowth2
FOREIGN KEY (pattern_class_edgegrowth2)
REFERENCES pattern_class_param (pattern_class);

ALTER TABLE pattern_class_output
ALTER COLUMN pattern_class_nogrowth DROP NOT NULL;



