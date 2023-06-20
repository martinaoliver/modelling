create table if not exists public.simulation_output(
    "simID" int not null,
    "model_param_id" text not null,
    "ssID" int not null ,
    "U_final_1D" numeric[],
    "U_final_2D" numeric[][],
    "U_record_1D" numeric[][],
    "U_record_2D" numeric[][][],

    primary key ("simID", "model_param_id", "ssID")
);


alter table simulation_output
    add constraint fk_simulation_output
        foreign key ("simID") references simulation_param ("simID");
--         foreign key (model_param)
