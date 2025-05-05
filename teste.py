from lib2to3.pgen2.pgen import DFAState
import requests
import time
from typing import List
from Bio import Entrez
import pandas as pd
import xml.etree.ElementTree as ET
import xmltodict
import time
import click
import json
import requests
import re
from urllib.parse import quote
from tqdm import tqdm

Entrez.email = 'pedro.rodrigues.rocha@usp.br'
Entrez.api_key = '94372c1e44aa8df940f2e7b790a761a25908'
DEFAULT_QUERY = 'metagenome AND aquatic'

def search_bioprojects(term, max_results):
    """Realiza pesquisa inicial de bioprojetos no NCBI"""
    print ("Iniciando a função search_bioprojects -(1)")
    print(f"Pesquisando bioprojetos com termo: {term}")
    handle = Entrez.esearch(db='bioproject', term=term, retmax=max_results)
    results = Entrez.read(handle)
    handle.close()
    print(f"Foram encontrados {results['Count']} resultados.")
    print("Fim da função search_bioprojects -(1)")
    return results['IdList']

def fetch_bioproject_data(id_list):
    """Obtém dados detalhados de cada bioprojeto"""
    print("Início da função fetch_bioproject_data -(2)")
    print(f"Obtendo dados para {len(id_list)} bioprojetos...")
    
    ids = []
    accessions = []
    organisms = []
    titles = []
    submissions = []
    
    for i, id_number in enumerate(id_list, 1):
        try:
            handle = Entrez.efetch(db='bioproject', id=id_number, rettype='xml', retmode='xml')
            print(f"Processando bioprojeto {i}/{len(id_list)} - ID: {id_number}")
            time.sleep(0.01)
            
            xml_data = handle.read().decode('utf-8')
            handle.close()
            
            data_dict = xmltodict.parse(xml_data)
            doc_summary = data_dict['RecordSet']['DocumentSummary']
            
            ids.append(doc_summary['@uid'])
            accessions.append(doc_summary['Project']['ProjectID']['ArchiveID']['@accession'])
            titles.append(doc_summary['Project']['ProjectDescr']['Title'])
            submissions.append(doc_summary['Submission']['@submitted'])
            
            try:
                organisms.append(doc_summary['Project']['ProjectType']['ProjectTypeSubmission']['Target']['Organism']['OrganismName'])
            except KeyError:
                organisms.append('NA')
                
        except Exception as e:
            print(f"Erro ao processar bioprojeto {id_number}: {e}")
            continue
    
    df = pd.DataFrame({
        'id': ids,
        'bioproject_id': accessions,
        'organism_name': organisms,
        'title': titles,
        'submission': submissions
    })
    
    print("\nDataFrame resultante da função fetch_bioproject_data:")
    print(df.head(25))  
    print(f"\nShape do DataFrame: {df.shape}")
    print("\nTipos de dados:")
    print(df.dtypes)
    print("Fim da função fetch_bioproject_data -(2)")
    return df

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import time
import random
from tqdm import tqdm
import pandas as pd

def setup_session():
    """Configure requests session with retry logic"""
    session = requests.Session()
    
    retry_strategy = Retry(
        total=5,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["HEAD", "GET", "OPTIONS"]
    )
    
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    
    return session

import pandas as pd
import time
import requests
from tqdm import tqdm

def link_bioproject_to_biosample(bioproject_ids):
    """
    Obtém relações completas BioProject ↔ BioSample com verificação bidirecional,
    e adiciona os accessions dos BioSamples ao dataframe parcial.
    
    Retorna:
    - DataFrame com todas as relações encontradas (df_relations)
    - Dicionário de relações inversas (biosample_to_projects)
    - DataFrame parcial com accessions dos BioSamples (df_parcial_relations)
    """
    session = requests.Session()
    parcial_relations = []
    relations = []
    biosample_to_projects = {}
    
    for bioproject_id in tqdm(bioproject_ids, desc="Processando BioProjects"):
        try:
            time.sleep(0.34)  # Respeita o limite de 3 requisições/segundo
            
            # Passo 1: Busca direta (BioProject → BioSample)
            response = session.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi",
                params={
                    "dbfrom": "bioproject",
                    "db": "biosample",
                    "id": bioproject_id,
                    "retmode": "json",
                    "email": "pedro.rodrigues.rocha@usp.br"
                },
                timeout=30
            )
            response.raise_for_status()
            data = response.json()
            
            for linkset in data["linksets"]:
                current_bioproject = linkset["ids"][0]
                for linksetdb in linkset.get("linksetdbs", []):
                    if linksetdb["dbto"] == "biosample":
                        for biosample_id in linksetdb["links"]:
                            # Extrai o accession (remove .X se existir)
                            biosample_accession = biosample_id.split('.')[0]
                            
                            # Adiciona relação direta
                            relations.append({
                                "BioProject_ID": current_bioproject,
                                "BioSample_ID": biosample_id,
                                "BioSample_Accession": biosample_accession
                            })
                            
                            # Adiciona ao dataframe parcial
                            parcial_relations.append({ 
                                "BioProject_ID": current_bioproject,
                                "BioSample_ID": biosample_id,
                                "BioSample_Accession": biosample_accession
                            })
                            
                            # Mapeamento inverso (BioSample → BioProjects)
                            if biosample_id not in biosample_to_projects:
                                biosample_to_projects[biosample_id] = set()
                            biosample_to_projects[biosample_id].add(current_bioproject)
                            
                            # Verificação inversa (BioSample → BioProject)
                            time.sleep(0.34)
                            try:
                                inverse_response = session.get(
                                    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi",
                                    params={
                                        "dbfrom": "biosample",
                                        "db": "bioproject",
                                        "id": biosample_id,
                                        "retmode": "json",
                                        "email": "pedro.rodrigues.rocha@usp.br"
                                    },
                                    timeout=30
                                )
                                inverse_response.raise_for_status()
                                inverse_data = inverse_response.json()
                                
                                for inv_linkset in inverse_data["linksets"]:
                                    for inv_linksetdb in inv_linkset.get("linksetdbs", []):
                                        if inv_linksetdb["dbto"] == "bioproject":
                                            for additional_project in inv_linksetdb["links"]:
                                                if additional_project != current_bioproject:
                                                    relations.append({
                                                        "BioProject_ID": additional_project,
                                                        "BioSample_ID": biosample_id,
                                                        "BioSample_Accession": biosample_accession
                                                    })
                                                    biosample_to_projects[biosample_id].add(additional_project)
                            
                            except Exception as inv_e:
                                print(f"Erro na verificação inversa para BioSample {biosample_id}: {str(inv_e)}")
                                continue
        
        except Exception as e:
            print(f"Erro no BioProject {bioproject_id}: {str(e)}")
            continue
    
    # Remove duplicatas e cria DataFrames
    df_relations = pd.DataFrame(relations).drop_duplicates()
    df_parcial_relations = pd.DataFrame(parcial_relations).drop_duplicates()
    
    # Salva em CSV 
    #df_relations.to_csv("0_bioproject_biosample_relations.csv", index=False)
    #df_parcial_relations.to_csv("0_bioproject_biosample_parcial.csv", index=False)

    # Obtendo a lista de biosamples
    lista_de_biosamples = df_parcial_relations['BioSample_ID'].tolist()
    
    return df_relations, biosample_to_projects, df_parcial_relations, lista_de_biosamples

from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd

def identifiers(lista_de_biosamples, email="pedro.rodrigues.rocha@usp.br"):
    """
    Busca identificadores de BioSample no NCBI para cada BioSample_ID.
    
    Parâmetros:
    lista_de_biosamples (list): Lista de IDs de BioSample
    email (str): E-mail para usar na API do NCBI
    
    Retorna:
    pandas.DataFrame: DataFrame com colunas 'BioSample_ID' e 'Accession'
    """
    Entrez.email = email
    results = []
    
    for biosample_id in lista_de_biosamples:
        try:
            # Obter os detalhes do BioSample
            handle = Entrez.efetch(db="biosample", id=biosample_id)
            biosample_data = handle.read().decode('utf-8')
            handle.close()
            
            # Parse do XML
            root = ET.fromstring(biosample_data)
            
            # Encontrar o elemento BioSample
            biosample = root.find('BioSample')
            accession = biosample.get('accession') if biosample is not None else None
            
            results.append({'BioSample_ID': biosample_id, 'Accession': accession})
                
        except Exception as e:
            print(f"Erro ao processar {biosample_id}: {str(e)}")
            results.append({'BioSample_ID': biosample_id, 'Accession': None})
    
    return pd.DataFrame(results)


def unindo_df_bioproject_biosample(df_bioproject, df_relations, df_parcial_relations, df_parcial_relations_atualizada, df_metadados):
    
    # Unindo as listas parciais
    parcial_relations = pd.merge(df_parcial_relations, df_parcial_relations_atualizada, on='BioSample_ID', how='left')
    parcial_relations.to_csv("0.5_parcial_relations.csv")

    # Unindo a relação bioproject -> biosample com os metadados de biosample
    df_metadados = df_metadados.rename(columns={"biosample_id": "Accession"})
    df_metadados_completos = pd.merge(parcial_relations, df_metadados, on='Accession', how='left')
    
    # Unindo a tabela completa de biosample com a tabela de bioproject
    df_bioproject = df_bioproject.rename(columns={"id": "BioProject_ID"})
    df_metadados_completos = pd.merge(df_bioproject, df_metadados_completos, left_on='BioProject_ID', right_on='BioProject_ID', how='right')
  
    # Unindo a tabela de bioproject/biosample com sra

    df_sra = df_sra.rename(columns={"BioProject_Origin": "BioProject_ID"})
    df_metadados_completos = pd.merge(df_metadados_completos, df_sra, left_on='BioProject_ID', right_on='BioProject_ID', how='right')

    return df_metadados_completos

def filter_data():
    """Filtra os dados conforme critérios especificados"""
    print("Iniciando a função filter_data -(6)")
    print("Filtrando dados...")
    
    # Filtros:

    
    print(f"Dados filtrados: {filtered.shape[0]} amostras restantes")
    print("A forma das primeiras linhas do df filtrado é:")
    print(filtered.head(25))
    print("Fim da função filter_data -(6)")
    return filtered

def save_and_summarize(df, filename):
    """Salva os dados e exibe um resumo"""
    print("Início da função save_and_summarize -(7)")
    print(f"Salvando dados em {filename}")
    df.to_csv(filename, index=False)
    
    print("\nResumo dos dados:")
    print("\nContagem por instrumento:")
    print(df['library_instrument'].value_counts())
    
    print("\nContagem por estratégia:")
    print(df['library_strategy'].value_counts())
    
    print("\nContagem por localização:")
    print(df['geo_loc_name'].value_counts())
    
    print("\nContagem por projeto:")
    print(df['bioproject_id'].value_counts())
    print("Fim da função save_and_summarize -(7)")
    
    #################################################################################
    #################################################################################
    #################################################################################
    
import pandas as pd
from Bio import Entrez
import xml.etree.ElementTree as ET
from time import sleep
from typing import Dict, Set, Union

Entrez.email = "pedro.rodrigues.rocha@usp.br"

def get_biosample_metadata(biosample_dict: Dict[str, Set[str]], batch_size: int = 100, delay: float = 1) -> pd.DataFrame:
    """
    Obtém metadados para BioSamples a partir de um dicionário de IDs
    
    Parâmetros:
    biosample_dict: Dicionário onde chaves são BioSample IDs e valores são conjuntos de IDs relacionados
    batch_size: Número de IDs para buscar por requisição
    delay: Atraso entre requisições (em segundos)
    
    Retorna:
    DataFrame: Metadados consolidados dos BioSamples principais
    """
    biosample_ids = list(biosample_dict.keys())
    all_data = []
    
    # Processar em lotes
    for i in range(0, len(biosample_ids), batch_size):
        batch = biosample_ids[i:i + batch_size]
        print(f"Processando lotes {i+1}-{i+len(batch)} de {len(biosample_ids)}...")
        
        try:
            # Buscar todos os IDs do lote de uma vez
            with Entrez.efetch(db="biosample", id=",".join(batch), retmode="xml") as handle:
                xml_data = handle.read()
            root = ET.fromstring(xml_data)
            
            for sample in root.findall('.//BioSample'):
                try:
                    attrs = sample.attrib
                    characteristics = {}
                    
                    for attr in sample.findall('.//Attribute'):
                        harmonized_name = attr.attrib.get('harmonized_name', '')
                        attribute_name = attr.attrib.get('attribute_name', harmonized_name)
                        characteristics[attribute_name] = attr.text
                    
                    # Adicionar os IDs relacionados ao registro
                    related_ids = biosample_dict.get(attrs.get('accession'), set())
                    sample_data = {
                        'biosample_id': attrs.get('accession'),
                        'related_ids': ", ".join(related_ids) if related_ids else None,
                        'publication_date': attrs.get('publication_date'),
                        'last_update': attrs.get('last_update'),
                        'submission_date': attrs.get('submission_date'),
                        **characteristics
                    }
                    all_data.append(sample_data)
                    
                except Exception as e:
                    print(f"Erro ao processar BioSample {attrs.get('accession', 'desconhecido')}: {str(e)}")
            
            sleep(delay)
            
        except Exception as e:
            print(f"Erro no lote {i//batch_size + 1}: {str(e)}")
            # Tentar processar os IDs individualmente
            for biosample_id in batch:
                try:
                    with Entrez.efetch(db="biosample", id=biosample_id, retmode="xml") as handle:
                        xml_data = handle.read()
                    root = ET.fromstring(xml_data)
                    sample = root.find('.//BioSample')
                    attrs = sample.attrib
                    characteristics = {}
                    
                    for attr in sample.findall('.//Attribute'):
                        harmonized_name = attr.attrib.get('harmonized_name', '')
                        attribute_name = attr.attrib.get('attribute_name', harmonized_name)
                        characteristics[attribute_name] = attr.text
                    
                    related_ids = biosample_dict.get(biosample_id, set())
                    all_data.append({
                        'biosample_id': biosample_id,
                        'related_ids': ", ".join(related_ids) if related_ids else None,
                        'publication_date': attrs.get('publication_date'),
                        'last_update': attrs.get('last_update'),
                        'submission_date': attrs.get('submission_date'),
                        **characteristics
                    })
                    sleep(delay/2)
                    
                except Exception as e:
                    print(f"Falha ao processar BioSample {biosample_id}: {str(e)}")
    
    # Criar DataFrame consolidado
    if all_data:
        df = pd.DataFrame(all_data)
        # Ordenar colunas (biosample_id primeiro)
        cols = ['biosample_id', 'related_ids', 'publication_date', 'last_update', 'submission_date'] + \
               [c for c in df.columns if c not in {'biosample_id', 'related_ids', 'publication_date', 'last_update', 'submission_date'}]
        return df[cols]
    
    return pd.DataFrame()

#########################################################################################################

import os
import time
import pandas as pd
from Bio import Entrez

def get_sra_metadata_for_multiple_bioprojects(bioproject_ids, email="pedro.rodrigues.rocha@usp.br"):
    """
    Obtém metadados SRA para uma lista de bioproject IDs e retorna:
    - DataFrame com metadados SRA (incluindo origem)
    - Lista de bioprojects sem registros SRA
    
    Args:
        bioproject_ids (list): Lista de IDs de bioprojects
        email (str): Email válido para uso na API do NCBI
    
    Returns:
        tuple: (pd.DataFrame com metadados, list de bioprojects sem SRA)
    """
    Entrez.email = email
    all_metadata = []
    bioprojects_without_sra = []
    delay = 1  # segundos entre requisições
    
    for i, bioproject_id in enumerate(bioproject_ids, 1):
        temp_file = None
        try:
            print(f"\nProcessando ({i}/{len(bioproject_ids)}): {bioproject_id}")
            
            # Buscar IDs SRA
            with Entrez.esearch(db="sra", term=f"{bioproject_id}[BioProject]", retmax=1000) as handle:
                record = Entrez.read(handle)
            
            sra_ids = record.get("IdList", [])
            
            if not sra_ids:
                print(f"  Nenhum SRA encontrado para {bioproject_id}")
                bioprojects_without_sra.append(bioproject_id)
                # Cria um DataFrame mínimo com apenas a coluna de origem
                empty_df = pd.DataFrame({'BioProject_Origin': [bioproject_id]})
                all_metadata.append(empty_df)
                continue
            
            print(f"  Encontrados {len(sra_ids)} registros SRA")
            
            # Obter metadados
            with Entrez.efetch(db="sra", id=",".join(sra_ids), rettype="runinfo", retmode="text") as handle:
                sra_metadata = handle.read().decode('utf-8')
            
            # Processar metadados
            temp_file = f"temp_{bioproject_id}.csv"
            with open(temp_file, "w", encoding='utf-8') as f:
                f.write(sra_metadata)
            
            df = pd.read_csv(temp_file)
            
            # Adicionar coluna de origem
            df.insert(0, 'BioProject_Origin', bioproject_id)
            
            all_metadata.append(df)
            
            time.sleep(delay)
            
        except Exception as e:
            print(f"  Erro no bioproject {bioproject_id}: {e}")
            # Registra como bioproject sem SRA e mantém o erro
            bioprojects_without_sra.append(bioproject_id)
            error_df = pd.DataFrame({
                'BioProject_Origin': [bioproject_id],
                'Error': [str(e)]
            })
            all_metadata.append(error_df)
            
        finally:
            # Garante que o arquivo temporário seja removido, mesmo se ocorrer erro
            if temp_file and os.path.exists(temp_file):
                os.remove(temp_file)
    
    # Consolidar resultados
    final_df = pd.concat(all_metadata, ignore_index=True) if all_metadata else pd.DataFrame()
    
    # Reorganizar colunas (se o DataFrame não estiver vazio)
    if not final_df.empty:
        cols = ['BioProject_Origin'] + [col for col in final_df.columns if col != 'BioProject_Origin']
        final_df = final_df[cols]
    final_df.to_csv("sra_metadata.csv")
    bioproject_sem_sra = pd.DataFrame(bioprojects_without_sra, columns=['BioProject_ID'])
    bioproject_sem_sra.to_csv("bioproject_sem_sra.csv")
    return final_df, bioprojects_without_sra
   
@click.command()
@click.option('--termo_da_pesquisa', default=DEFAULT_QUERY, help="Termo de pesquisa no NCBI.")
@click.option('--num_resultados', default=100, type=int, help="Número máximo de resultados.")
@click.option('--termo_estrategia', default="wgs", help="Tipo de estratégia de sequenciamento.")
@click.option('--termo_instrumento', default="metagenomic,genomic", help="Tipos de instrumentos aceitos.")
@click.option('--termo_localizacao', default="brazil", help="Localização geográfica das amostras.")
def main(termo_da_pesquisa, num_resultados, termo_estrategia, termo_instrumento, termo_localizacao):
    """Pipeline principal para busca e filtragem de dados do NCBI"""
    
    # Processa parâmetros
    instrumentos = [x.strip() for x in termo_instrumento.split(',')]
    localizacoes = [x.strip() for x in termo_localizacao.split(',')]
    
    # 1. Busca bioprojetos
    bioproject_ids = search_bioprojects(termo_da_pesquisa, num_resultados)
    print(bioproject_ids)
    
    # 2. Obtém dados dos bioprojetos
    df_bioproject = fetch_bioproject_data(bioproject_ids)
    df_bioproject.to_csv("1_df_bioproject.csv")
    
    # 3. Encontra biosamples relacionados
    df_relations, biosample_ids, df_parcial_relations, lista_de_biosamples = link_bioproject_to_biosample(bioproject_ids)

    # 3.5. Identificadores
    df_parcial_relations_atualizada = identifiers(lista_de_biosamples, email="pedro.rodrigues.rocha@usp.br")
    
    # 4. Faz a pesquisa dos metadados dos biosamples
    df_metadados = get_biosample_metadata(biosample_ids, batch_size=10)
    df_metadados.to_csv("1_df_metadados.csv")
    # 5. Obtenção da metadata dos sra
    df_sra, bioprojects_sem_sra = get_sra_metadata_for_multiple_bioprojects(bioproject_ids, email="pedro.rodrigues.rocha@usp.br")
    print(df_sra)
    
    # 5. Combina dados de bioprojeto e biosample
    df_metadados_completos = unindo_df_bioproject_biosample(df_bioproject, df_relations, df_parcial_relations, df_parcial_relations_atualizada, df_metadados)
    
    # 6. Filtra os dados
    #filtered_df = filter_data(final_df, instrumentos, [termo_estrategia], localizacoes)
    
    # 7. Salva e exibe resumo
    #save_and_summarize(filtered_df, 'lista_ncbi_filtrada.csv')

    df_metadados_completos.to_csv("finalmente.csv")

if __name__ == "__main__":
    main()
