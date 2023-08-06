create table if not exists public.pattern_class_param(
    "pattern_class" text not null,
    "converged" boolean,
    "flat" boolean,
    "regular" boolean,
    "peak_number" int,
    primary key ("pattern_class")
);

alter table
INSERT INTO pattern_class_param ("pattern_class", "flat","converged")
VALUES ('homogeneous', True, True);
INSERT INTO pattern_class_param ("pattern_class", "flat","converged")
VALUES ('temporal oscillator', True, False);
INSERT INTO pattern_class_param ("pattern_class", "flat","converged", "regular")
VALUES ('stationary regular pattern', False, True, True);
INSERT INTO pattern_class_param ("pattern_class", "flat","converged", "regular")
VALUES ('stationary irregular pattern', False, True, False);

alter table pattern_class_param rename column peak_number to n_peaks;

INSERT INTO pattern_class_param("pattern_class", "n_peaks")
VALUES ('no pattern, homogeneous',1 );
INSERT INTO pattern_class_param("pattern_class", "n_peaks")
VALUES ('no pattern, boundary effect' , 2 );
INSERT INTO pattern_class_param("pattern_class", "n_peaks")
VALUES ('weak pattern' , 3 );
INSERT INTO pattern_class_param("pattern_class", "n_peaks")
VALUES ('intermediate pattern' , 4 );
INSERT INTO pattern_class_param("pattern_class", "n_peaks")
VALUES ('strong pattern' , 5 );

alter table pattern_class_param drop column regular;

