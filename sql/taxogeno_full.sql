-- CREATE EXTENSION IF NOT EXISTS pgcrypto;

-----------+----------+
-- ESQUEMA | taxogeno |
---------- +----------+
-- Este esquema es para agrupar las tablas que están medio normalizadas.
-- Dentro de las tablas que contiene hay que meter la información normalizada,
-- sacándola de donde se pueda. En el fichero DML lo que se hace es una select
-- sobre la tabla bruta bulkinput, para hacer un insert masivo con los datos
-- que interesen.
-----------------------
drop schema if exists taxogeno cascade;
	      create schema taxogeno;

----------------------------------
--  esquema    tabla
--+----------+-------------------+
--| taxogeno | taxogeno.taxonomy |
--+----------+-------------------+
CREATE TABLE taxogeno.taxonomy(
  ncbitaxid        integer PRIMARY KEY,
  parent_ncbitaxid integer not null,
  node_rank        character varying,
  scientific_name  character varying not null,
  is_ancestor      boolean not null
);

-----------------------------------
--  esquema    tabla
--+----------+--------------------+
--| taxogeno | taxogeno.proteome  |
--+----------+--------------------+
CREATE TABLE taxogeno.proteome(
  proteomeid serial PRIMARY KEY,

  -- -- -- -- -- -- -- -- -- -- -- -- -- --
  basename         character varying,
  is_userproteome  boolean   not null default FALSE,
  gcaid		   character varying,
  ncbitaxid        integer references taxogeno.taxonomy(ncbitaxid),
  dbname           character varying,
  do_updatencbiassemblyinfo boolean default TRUE not null,
  sourceurl        character varying,

  -- -- -- -- -- -- -- -- -- -- -- -- -- --
  multifastamd5sum        uuid,
  UNIQUE (multifastamd5sum, is_userproteome),
  sma3sannotationmd5sum   uuid not null,
  UNIQUE (sma3sannotationmd5sum, is_userproteome),

  -- -- -- -- -- -- -- -- -- -- -- -- -- --
  creationtimestamp       timestamp default now(),

  genecount		  integer
  -- genecountaaseq	  integer,
  -- genecountgocannot	  integer,
  -- genecountgofannot	  integer,
  -- genecountgopannot	  integer,
  -- genecountgoslimannot	  integer,
  -- genecountkeywordannot	  integer,
  -- genecountenzymeannot	  integer,
  -- genecountpathwayannot	  integer
);

create table taxogeno.proteome_tag_rel(
  proteomeid integer references taxogeno.proteome(proteomeid),
  tag character varying,
  primary key (proteomeid,tag)
);

-------------------------------
--  esquema    tabla
--+----------+----------------+
--| taxogeno | taxogeno.gene  |
--+----------+----------------+
create table taxogeno.gene(
  geneid          serial        primary key,
  
  -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  proteomeid       integer       references taxogeno.proteome(proteomeid) not null,
  fastaheader      character varying not null,
  unique (proteomeid,fastaheader),
  fastashortheader character varying not null,	
  unique (proteomeid,fastashortheader),

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  genename         character varying,
  genedescription  text,
  aasequence       text
);
CREATE INDEX gene_proteomeid_rel_idx ON taxogeno.gene USING btree (proteomeid);

-- GO P
--------------------------------------
--  esquema    tabla
--+----------+-----------------------+
--| taxogeno | taxogeno.gene_gop_rel |
--+----------+-----------------------+
-- Esta tabla es como si fuera una tabla de relación pero sin serlo.
-- No apunta con foreign keys a ninguna tabla que contenga los id de GeneOntology, aunque existe.
-- Los GO se añaden, se deprecan...
-- Controlar su integridad no es tarea de esta base de datos. -> owltools
-- Política de no hacer foreign key a bases de datos externas.
--------------------------------------
create table taxogeno.gene_gop_rel(
  geneid integer references taxogeno.gene(geneid),
  goid char(10),
  primary key(geneid,goid)
);

-- GO F
--------------------------------------
--  esquema    tabla
--+----------+-----------------------+
--| taxogeno | taxogeno.gene_gof_rel |
--+----------+-----------------------+
-- Esta tabla es como si fuera una tabla de relación pero sin serlo.
-- No apunta con foreign keys a ninguna tabla que contenga los id de GeneOntology, aunque existe.
-- Los GO se añaden, se deprecan...
-- Controlar su integridad no es tarea de esta base de datos. -> owltools
-- Política de no hacer foreign key a bases de datos externas.
--------------------------------------
create table taxogeno.gene_gof_rel(
  geneid integer references taxogeno.gene(geneid),
  goid char(10),
  primary key(geneid,goid)
);

-- GO C
--------------------------------------
--  esquema    tabla
--+----------+-----------------------+
--| taxogeno | taxogeno.gene_goc_rel |
--+----------+-----------------------+
-- Esta tabla es como si fuera una tabla de relación pero sin serlo.
-- No apunta con foreign keys a ninguna tabla que contenga los id de GeneOntology, aunque existe.
-- Los GO se añaden, se deprecan...
-- Controlar su integridad no es tarea de esta base de datos. -> owltools
-- Política de no hacer foreign key a bases de datos externas.
--------------------------------------
create table taxogeno.gene_goc_rel(
  geneid integer references taxogeno.gene(geneid),
  goid char(10),
  primary key(geneid,goid)
);

-- GOSLIM
--------------------------------------
--  esquema    tabla
--+----------+-----------------------+
--| taxogeno | taxogeno.gene_goslim_rel |
--+----------+-----------------------+
-- Esta tabla es como si fuera una tabla de relación pero sin serlo.
-- No apunta con foreign keys a ninguna tabla que contenga los id de GeneOntology, aunque existe.
-- Los GO se añaden, se deprecan...
-- Controlar su integridad no es tarea de esta base de datos. -> owltools
--------------------------------------
create table taxogeno.gene_goslim_rel(
  geneid integer references taxogeno.gene(geneid),
  goid char(10),
  primary key(geneid,goid)
);

----------
-- taxogeno.gene_keyword_rel
----------
create table taxogeno.gene_keyword_rel(
  geneid integer references taxogeno.gene(geneid),
  keyword varchar(100),
  primary key(geneid,keyword)
);

----------
-- taxogeno.gene_enzyme_rel
----------
-- Política de no hacer foreign key a bases de datos externas
create table taxogeno.gene_enzyme_rel(
  geneid integer references taxogeno.gene(geneid),
  ec varchar(100),
  primary key(geneid,ec)
);

----------
-- taxogeno.gene_pathway_rel
----------
-- Política de no hacer foreign key a bases de datos externas
create table taxogeno.gene_pathway_rel(
  geneid integer references taxogeno.gene(geneid),
  pathway varchar(100),
  primary key(geneid,pathway)
);

-- SUMMARIES
-- Generados por owltools
create table taxogeno.generated_goslim_summary(
  proteomeid integer references taxogeno.proteome(proteomeid),
  annotkwid char(10),
  annotkwcount integer not null,
  primary key (proteomeid,annotkwid)
);
create table taxogeno.generated_keyword_summary(
  proteomeid integer references taxogeno.proteome(proteomeid),
  annotkwid char(10),
  annotkwcount integer not null,
  primary key (proteomeid,annotkwid)
);
create table taxogeno.keyword_summary(
  proteomeid integer references taxogeno.proteome(proteomeid),
  annotkwid char(10),
  annotkwcount integer not null,
  primary key (proteomeid,annotkwid)
);
-- Índices extra para la tabla taxogeno.generated_summary
----------------------------------------------------
-- Hay campos que deberían tener un índice pues van a hacerse búsquedas y joins
-- con igualdades en ellos.
------

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
CREATE TABLE taxogeno.euclidean_distances_generated_goslim (
       proteomeid_greatest integer references taxogeno.proteome(proteomeid),
       proteomeid_least integer references taxogeno.proteome(proteomeid),
       distance numeric not null,
       check (proteomeid_greatest>proteomeid_least),
       primary key (proteomeid_greatest,proteomeid_least)
);
CREATE TABLE taxogeno.euclidean_distances_keyword (
       proteomeid_greatest integer references taxogeno.proteome(proteomeid),
       proteomeid_least integer references taxogeno.proteome(proteomeid),
       distance numeric not null,
       check (proteomeid_greatest>proteomeid_least),
       primary key (proteomeid_greatest,proteomeid_least)
);


------------------
------------------


CREATE OR REPLACE PROCEDURE taxogeno.delete_proteome(p_proteomeid integer) AS $$
  DECLARE
  v_geneid integer;
BEGIN
  DELETE FROM taxogeno.generated_goslim_summary
   where proteomeid=p_proteomeid;
  DELETE FROM taxogeno.euclidean_distances_keyword
   where proteomeid_greatest=p_proteomeid OR proteomeid_least=p_proteomeid;
  DELETE FROM taxogeno.euclidean_distances_generated_goslim
   where proteomeid_greatest=p_proteomeid OR proteomeid_least=p_proteomeid;
  FOR v_geneid in (SELECT geneid FROM taxogeno.gene WHERE proteomeid=p_proteomeid) LOOP
    DELETE FROM taxogeno.gene_goc_rel WHERE geneid=v_geneid;
    DELETE FROM taxogeno.gene_gof_rel WHERE geneid=v_geneid;
    DELETE FROM taxogeno.gene_gop_rel WHERE geneid=v_geneid;
    DELETE FROM taxogeno.gene_goslim_rel WHERE geneid=v_geneid;
    DELETE FROM taxogeno.gene_keyword_rel WHERE geneid=v_geneid;
    DELETE FROM taxogeno.gene_pathway_rel WHERE geneid=v_geneid;
    DELETE FROM taxogeno.gene_enzyme_rel WHERE geneid=v_geneid;
    DELETE FROM taxogeno.gene WHERE geneid=v_geneid;
  END LOOP;
  
  DELETE FROM taxogeno.proteome WHERE proteomeid=p_proteomeid;  
END;
$$ LANGUAGE PLPGSQL;

CREATE OR REPLACE PROCEDURE taxogeno.delete_all_proteomes() AS $$
  DECLARE
  v_geneid integer;
BEGIN
  DELETE FROM taxogeno.generated_goslim_summary;
  DELETE FROM taxogeno.euclidean_distances_generated_goslim;
  DELETE FROM taxogeno.gene_goc_rel;
  DELETE FROM taxogeno.gene_gof_rel;
  DELETE FROM taxogeno.gene_gop_rel;
  DELETE FROM taxogeno.gene_goslim_rel;
  DELETE FROM taxogeno.gene_keyword_rel;
  DELETE FROM taxogeno.gene_pathway_rel;
  DELETE FROM taxogeno.gene_enzyme_rel;
  DELETE FROM taxogeno.gene;
  DELETE FROM taxogeno.proteome;  
END;
$$ LANGUAGE PLPGSQL;


  CREATE OR REPLACE FUNCTION taxogeno.taxonomy_jsonb_children(p_ncbitaxid integer)
    RETURNS jsonb LANGUAGE PLPGSQL AS $$
    DECLARE
    v_jsonbobj jsonb;
    v_ncbitaxid integer;
  BEGIN
    SELECT 
      jsonb_build_object('ncbitaxid', tax1.ncbitaxid,
			 'parent_ncbitaxid', tax1.parent_ncbitaxid,
			 'node_rank', tax1.node_rank,
			 'scientific_name', tax1.scientific_name,
			 'children',(SELECT jsonb_agg(taxogeno.taxonomy_jsonb_children(ncbitaxid))
				       FROM taxogeno.taxonomy tax2
				      WHERE tax2.parent_ncbitaxid=tax1.ncbitaxid
					AND tax2.parent_ncbitaxid<>tax2.ncbitaxid))		   
      INTO v_jsonbobj
      FROM taxogeno.taxonomy tax1
     WHERE ncbitaxid=p_ncbitaxid;
    RETURN v_jsonbobj;
  END;
  $$;

CREATE OR REPLACE PROCEDURE taxogeno.insert_taxon(p_ncbitaxid integer)     
AS
$$
  DECLARE
    v_is_ancestor boolean;
    v_ncbitaxid integer;
    v_parent_ncbitaxid integer;
    v_node_rank varchar(100);
    v_scientific_name varchar(100);
    v_pretaxonomy_rec record;

    v_inserted_taxons integer;
  BEGIN
  
  FOR v_pretaxonomy_rec in (SELECT * FROM biosql.taxon_ancestors(p_ncbitaxid) ta)
  LOOP
    --BEGIN VARIABLE EXPANSION
    IF v_pretaxonomy_rec.ncbi_taxon_id<>p_ncbitaxid THEN
      v_is_ancestor   := TRUE;
    ELSE
      v_is_ancestor   := FALSE;
    END IF;
    v_ncbitaxid       := v_pretaxonomy_rec.ncbi_taxon_id;
    v_node_rank       := v_pretaxonomy_rec.node_rank;
    v_scientific_name := v_pretaxonomy_rec.name;
    -- END VARIABLE EXPANSION

    -- BEGIN PARENT SELECTION
    SELECT ncbi_taxon_id INTO v_parent_ncbitaxid
     FROM biosql.taxon
    WHERE biosql.taxon.taxon_id = v_pretaxonomy_rec.parent_taxon_id;
    -- END PARENT SELECTION

    -- BEGIN BUFFER TABLE INSERTION
    INSERT INTO taxogeno.taxonomy(ncbitaxid,   parent_ncbitaxid,   node_rank,   scientific_name,   is_ancestor)
                          VALUES (v_ncbitaxid, v_parent_ncbitaxid, v_node_rank, v_scientific_name, v_is_ancestor)
    ON CONFLICT (ncbitaxid)
    DO UPDATE
    SET parent_ncbitaxid=EXCLUDED.parent_ncbitaxid,
        node_rank=EXCLUDED.node_rank,
	scientific_name=EXCLUDED.scientific_name,
	is_ancestor=EXCLUDED.is_ancestor;
    -- END BUFFER TABLE INSERTION
  END LOOP;

  -- SELECT count(1) into v_inserted_taxons from taxogeno.taxon_ancestors(p_ncbitaxid);
  -- IF v_inserted_taxons<1 THEN
  --   RAISE EXCEPTION USING
  --   ERRCODE='P0002';
  -- END IF;

END;
$$ LANGUAGE plpgsql;


CREATE OR REPLACE PROCEDURE taxogeno.delete_taxon(p_ncbitaxid integer)
AS
$$
  DECLARE
    v_remaining_ncbitaxid integer;
    v_remaining_ncbitaxid_array_step integer[];
    v_remaining_ncbitaxid_array integer[];
    v_target_ncbitaxid_array integer[];
  BEGIN
    select array_agg(ncbitaxid) into v_target_ncbitaxid_array from taxogeno.taxon_ancestors(p_ncbitaxid);
    
    FOR v_remaining_ncbitaxid IN (SELECT ncbitaxid
				    FROM taxogeno.taxonomy
				   WHERE is_ancestor IS FALSE
				     AND ncbitaxid <> p_ncbitaxid) LOOP
      
      select array_agg(ncbitaxid)
	into v_remaining_ncbitaxid_array_step
	from taxogeno.taxon_ancestors(v_remaining_ncbitaxid);

      v_remaining_ncbitaxid_array:= v_remaining_ncbitaxid_array || v_remaining_ncbitaxid_array_step;
	       
    END LOOP;

    delete from taxogeno.taxonomy
     where ncbitaxid <> any(v_remaining_ncbitaxid_array)
       and ncbitaxid  = any(v_target_ncbitaxid_array);

  END;
$$ LANGUAGE plpgsql;


CREATE OR REPLACE FUNCTION taxogeno.taxon_ancestors(integer)
RETURNS TABLE(ncbitaxid integer, parent_ncbitaxid integer, node_rank character varying, scientific_name character varying, is_ancestor boolean)
LANGUAGE sql
AS $function$
	WITH RECURSIVE
	--  ancestors
	rec_a (ncbitaxid, parent_ncbitaxid, node_rank, scientific_name, is_ancestor) AS
	(
		SELECT
	            taxogeno.taxonomy.ncbitaxid,
		    taxogeno.taxonomy.parent_ncbitaxid,
		    taxogeno.taxonomy.node_rank, taxogeno.taxonomy.scientific_name,
		    taxogeno.taxonomy.is_ancestor
		 FROM taxogeno.taxonomy
		WHERE taxogeno.taxonomy.ncbitaxid=$1 --OBJETIVO
		
		UNION ALL
		
		SELECT
		   taxogeno.taxonomy.ncbitaxid,
		   taxogeno.taxonomy.parent_ncbitaxid,
		   taxogeno.taxonomy.node_rank,
		   taxogeno.taxonomy.scientific_name,
		   taxogeno.taxonomy.is_ancestor
		 FROM rec_a, taxogeno.taxonomy
		WHERE taxogeno.taxonomy.ncbitaxid = rec_a.parent_ncbitaxid
		  AND rec_a.ncbitaxid<>rec_a.parent_ncbitaxid
	)
	SELECT ncbitaxid, parent_ncbitaxid, node_rank, scientific_name, is_ancestor FROM rec_a;
$function$;

CREATE OR REPLACE FUNCTION taxogeno.taxon_descendants(integer)
RETURNS TABLE(ncbitaxid integer, parent_ncbitaxid integer, node_rank character varying, scientific_name character varying, is_ancestor boolean)
LANGUAGE sql                                             
AS $function$
	WITH RECURSIVE
	-- descendants
	rec_d (ncbitaxid, parent_ncbitaxid, node_rank, scientific_name, is_ancestor) AS
	(
		SELECT
		    taxogeno.taxonomy.ncbitaxid,
		    taxogeno.taxonomy.parent_ncbitaxid,
		    taxogeno.taxonomy.node_rank,
		    taxogeno.taxonomy.scientific_name,
		    taxogeno.taxonomy.is_ancestor
		 FROM taxogeno.taxonomy
		WHERE taxogeno.taxonomy.ncbitaxid=$1 --OBJETIVO
		
		UNION ALL
		
		SELECT
		    taxogeno.taxonomy.ncbitaxid,
		    taxogeno.taxonomy.parent_ncbitaxid,
		    taxogeno.taxonomy.node_rank,
		    taxogeno.taxonomy.scientific_name,
		    taxogeno.taxonomy.is_ancestor
		 FROM rec_d, taxogeno.taxonomy
		WHERE taxogeno.taxonomy.parent_ncbitaxid = rec_d.ncbitaxid
		  AND rec_d.ncbitaxid<>rec_d.parent_ncbitaxid
	)
	SELECT ncbitaxid, parent_ncbitaxid, node_rank, scientific_name, is_ancestor FROM rec_d;
$function$;
-- END procedures taxogeno taxonomy

-------------------------------
-- UNIPROT
-------------------------------
drop schema if exists uniprot cascade;
create schema uniprot;

--------------------------------
-- Keyword catalog
---------------------------------
CREATE TABLE uniprot.uniprot_keyword_category(
  keywordcategory     character varying primary key
);
CREATE TABLE uniprot.uniprot_keyword(
  keywordid           character varying primary key, -- KWXXXXXXX
  keyword             character varying unique not null, 
  keyworddescription  character varying,
  keywordcategory     character varying not null references uniprot.uniprot_keyword_category
);

--------------------------------
-- Proteomes catalog
---------------------------------
create table uniprot.proteomes_catalog(
  upid character varying,
  gcaid character varying,
  ncbitaxid integer,
  busco character varying,
  cpd character varying
);

--------------------------------
-- Pathway catalog
---------------------------------
-- [ ] por hacer


-------------------------------
-- NCBI GenBank
-------------------------------
drop schema if exists ncbi cascade;
create schema ncbi;
create table ncbi.assembly_summary_genbank(
  assembly_accession character varying,
  bioproject character varying,
  biosample character varying,
  wgs_master character varying,
  refseq_category character varying,
  taxid integer,
  species_taxid integer,
  organism_name character varying,
  infraspecific_name character varying,
  isolate character varying,
  version_status character varying,
  assembly_level character varying,
  release_type character varying,
  genome_rep character varying,
  seq_rel_date character varying,
  asm_name character varying,
  submitter character varying,
  gbrs_paired_asm character varying,
  paired_asm_comp character varying,
  ftp_path character varying,
  excluded_from_refseq character varying,
  relation_to_type_material character varying
);



DROP SCHEMA IF EXISTS biosql cascade;
CREATE SCHEMA biosql;

CREATE OR REPLACE FUNCTION biosql.taxon_ancestors(integer)
 RETURNS TABLE(ncbi_taxon_id integer, node_rank character varying, taxon_id integer, parent_taxon_id integer, name text)
 LANGUAGE sql
AS $function$
	WITH RECURSIVE
	--  ancestors
	rec_a (ncbi_taxon_id, node_rank, taxon_id, parent_taxon_id, name) AS
	(
		SELECT biosql.taxon.ncbi_taxon_id, biosql.taxon.node_rank, biosql.taxon.taxon_id, biosql.taxon.parent_taxon_id , biosql.taxon_name.name
		FROM biosql.taxon
		inner join biosql.taxon_name
		on (biosql.taxon_name.taxon_id= biosql.taxon.taxon_id and biosql.taxon_name.name_class='scientific name')
		WHERE biosql.taxon.ncbi_taxon_id=$1 --OBJETIVO
		UNION ALL
		SELECT biosql.taxon.ncbi_taxon_id, biosql.taxon.node_rank, biosql.taxon.taxon_id, biosql.taxon.parent_taxon_id , biosql.taxon_name.name
		FROM rec_a, biosql.taxon
		inner join biosql.taxon_name
		on (biosql.taxon_name.taxon_id= biosql.taxon.taxon_id  and biosql.taxon_name.name_class='scientific name')
		WHERE biosql.taxon.taxon_id = rec_a.parent_taxon_id
		AND rec_a.ncbi_taxon_id<>1
	)
	SELECT ncbi_taxon_id, node_rank, taxon_id,parent_taxon_id,name FROM rec_a;
$function$;

CREATE OR REPLACE FUNCTION biosql.taxon_descendants(integer)
RETURNS TABLE(ncbi_taxon_id integer, node_rank varchar(50), taxon_id integer, parent_taxon_id integer, name text)
LANGUAGE sql                                           
AS $function$
	WITH RECURSIVE
	-- descendants
	rec_d (ncbi_taxon_id, node_rank, taxon_id,parent_taxon_id, name) AS
	(
		SELECT taxon.ncbi_taxon_id, taxon.node_rank, taxon.taxon_id, taxon.parent_taxon_id , taxon_name.name
		FROM biosql.taxon
		inner join biosql.taxon_name
		on (taxon_name.taxon_id= taxon.taxon_id and taxon_name.name_class='scientific name')
		WHERE taxon.ncbi_taxon_id=$1 --OBJETIVO
		UNION ALL
		SELECT taxon.ncbi_taxon_id, taxon.node_rank, taxon.taxon_id, taxon.parent_taxon_id , taxon_name.name
		FROM rec_d, biosql.taxon
		inner join biosql.taxon_name
		on (taxon_name.taxon_id=taxon.taxon_id and taxon_name.name_class='scientific name')
		where taxon.parent_taxon_id = rec_d.taxon_id
		AND rec_d.ncbi_taxon_id<>1
	)
	SELECT ncbi_taxon_id, node_rank,taxon_id,parent_taxon_id,name FROM rec_d;
$function$;
  
CREATE OR REPLACE FUNCTION biosql.get_ncbitaxid_by_name(p_name character varying)
RETURNS integer 
AS $$
	DECLARE
		v_ncbitaxid integer;
	BEGIN
		SELECT DISTINCT ncbi_taxon_id INTO v_ncbitaxid
		FROM biosql.taxon_name
		INNER JOIN biosql.taxon ON (biosql.taxon_name.taxon_id=biosql.taxon.taxon_id)
		WHERE biosql.taxon_name.name=p_name;

		RETURN v_ncbitaxid;
	END;
  $$ LANGUAGE plpgsql;

SET search_path to public;
  
-------------------------------
-- GENE ONTOLOGY
-------------------------------
drop schema if exists gene_ontology CASCADE;
create schema gene_ontology;
create table gene_ontology.gene_ontology_xml(filename character varying primary key, xmldata xml);

create table gene_ontology.gene_ontology (
  goid     char(10) primary key,
  golabel  character varying not null,
  goaspect character varying
);

create table gene_ontology.goslim_generic (
  goid     char(10) references gene_ontology.gene_ontology(goid) primary key,
  golabel  character varying not null,
  goaspect character varying
);
