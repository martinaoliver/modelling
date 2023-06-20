create table if not exists public.simulation_param(
    "simID" int not null,
    "L" numeric,
    "T" numeric,
    "J" numeric,
    "N" numeric,
    "dx" numeric,
    "dt" numeric,
    "boundaryCoeff" numeric,
    "growth" text,
    "shape" text,
    "growth rate" numeric,
    primary key ("simID")
);




