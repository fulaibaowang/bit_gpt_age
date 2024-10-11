# ---
# jupyter:
#   jupytext:
#     cell_markers: '{{{,}}}'
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.1
#   kernelspec:
#     display_name: Python 3.10.8
#     language: python
#     name: py3.10.8
# ---

from llama_index import VectorStoreIndex, SimpleDirectoryReader

from llama_index import StorageContext, load_index_from_storage

storage_context = StorageContext.from_defaults(persist_dir="/nexus/posix0/MAGE-flaski/service/projects/data/Bioinformatics/bit_gpt_age/pmc_papers/")
index = load_index_from_storage(storage_context)

help(load_index_from_storage)

help(VectorStoreIndex)

help(StorageContext.from_defaults)

help(StorageContext)

from llama_index import VectorStoreIndex, SimpleDirectoryReader, StorageContext, load_index_from_storage

help(VectorStoreIndex.from_documents)

import llama_index
help( llama_index.storage.storage_context.StorageContext )

help(llama_index.indices.service_context.ServiceContext)

embed_model = LangchainEmbedding(HuggingFaceEmbeddings(model_name='TheBloke/Llama-2-13B-chat-GGML'))

# {{{
from llama_index import LLMPredictor, ServiceContext
from langchain.embeddings.huggingface import HuggingFaceEmbeddings
from llama_index import LangchainEmbedding, ServiceContext
from llama_index import VectorStoreIndex, SimpleDirectoryReader, StorageContext, load_index_from_storage


embed_model = LangchainEmbedding(HuggingFaceEmbeddings(model_name='BAAI/bge-small-en'))
service_context = ServiceContext.from_defaults(embed_model=embed_model)
documents = SimpleDirectoryReader('/flaski_private/pmc_papers', recursive=True ).load_data()
index = VectorStoreIndex.from_documents(documents, service_context=service_context)

# llm_predictor = LLMPredictor(llm=ChatOpenAI(temperature=0, model_name="gpt-3.5-turbo"))

# }}}

help(StorageContext)

help(ServiceContext.from_defaults)

from langchain import OpenAI
llm_predictor = LLMPredictor(llm=OpenAI(temperature=0, model_name="llama2-13b-chat"))

help(LangchainEmbedding)

help(HuggingFaceEmbeddings)

help(VectorStoreIndex.from_documents)

help(VectorStoreIndex)


