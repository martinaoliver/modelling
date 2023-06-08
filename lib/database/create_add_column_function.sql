
-- use like: SELECT f_add_col('public.kat', 'my_column_name', 'sql_data_type');
CREATE OR REPLACE function f_add_col(_tbl regclass, _col  text, _type regtype)
  RETURNS bool
  LANGUAGE plpgsql AS
$func$
BEGIN
   IF EXISTS (SELECT FROM pg_attribute
              WHERE  attrelid = _tbl
              AND    attname  = _col
              AND    NOT attisdropped) THEN
      RETURN false;
   ELSE
      EXECUTE format('ALTER TABLE %s ADD COLUMN "%I" %s', _tbl, _col, _type);
      RETURN true;
   END IF;
END
$func$;