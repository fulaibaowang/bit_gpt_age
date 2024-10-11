# bit_gpt_age


## packages

```
pip3 install llama-index llama-cpp-python=0.1.78 sentence_transformers --user

```

## traninging own data

```
from llama_index import VectorStoreIndex, SimpleDirectoryReader, StorageContext, load_index_from_storage
documents = SimpleDirectoryReader('/nexus/posix0/MAGE-flaski/service/projects/data/Bioinformatics/bit_gpt_age/pmc_papers', recursive=True ).load_data()
index = VectorStoreIndex.from_documents(documents, show_progress=True)
```

## storing 
```
index.storage_context.persist("/nexus/posix0/MAGE-flaski/service/projects/data/Bioinformatics/bit_gpt_age/trained/") 
```

## loading from query data
```
storage_context = StorageContext.from_defaults(persist_dir="/nexus/posix0/MAGE-flaski/service/projects/data/Bioinformatics/bit_gpt_age/trained/")
index = load_index_from_storage(storage_context)
```

## querying
```
query_engine = index.as_query_engine()
response = query_engine.query("What is a mitochondria?")
print(response) # or response.response
```

## papers

using `pubid_parser.py` and PR/PURE provided list of papers to download publicly available papers.

# Load in a specific embedding model

**https://blog.futuresmart.ai/mastering-llamaindex-create-save-load-indexes-customize-llms-prompts-embeddings**


```
from llama_index import VectorStoreIndex, SimpleDirectoryReader, StorageContext, load_index_from_storage
from langchain.embeddings.huggingface import HuggingFaceEmbeddings
from llama_index import LangchainEmbedding, ServiceContext

embed_model = LangchainEmbedding(HuggingFaceEmbeddings(model_name='BAAI/bge-small-en')) # at this step the embedding files are downloaded
service_context = ServiceContext.from_defaults(embed_model=embed_model)
documents = SimpleDirectoryReader('/flaski_private/pmc_papers', recursive=True ).load_data()
index = VectorStoreIndex.from_documents(documents, service_context=service_context)
```

flaski input: /flaski_private/pmc_papers

flaski output: /flaski_private/agebot/





index.storage_context.persist("/flaski_private/agebot/") 

pip3 install llama-index==0.8.23.post1 llama-cpp-python==0.1.78 sentence_transformers==2.2.2 --user

from llama_index import VectorStoreIndex, SimpleDirectoryReader, StorageContext, load_index_from_storage
from langchain.embeddings.huggingface import HuggingFaceEmbeddings
from llama_index import LangchainEmbedding, ServiceContext
storage_context = StorageContext.from_defaults(persist_dir="/flaski_private/agebot/")
index = load_index_from_storage(storage_context)


https://huggingface.co/TheBloke/Llama-2-13B-chat-GGML/resolve/main/llama-2-13b-chat.ggmlv3.q4_0.bin to path /tmp/llama_index/models/llama-2-13b-chat.ggmlv3.q4_0.bin


mkdir -p /tmp/llama_index/models/ && cp /flaski_private/llama-2-13b-chat.ggmlv3.q4_0.bin /tmp/llama_index/models/
pip3 install llama-index llama-cpp-python==0.1.78 sentence_transformers --user
python3
from llama_index import VectorStoreIndex, SimpleDirectoryReader, StorageContext, load_index_from_storage
from langchain.embeddings.huggingface import HuggingFaceEmbeddings
from llama_index import LangchainEmbedding, ServiceContext
embed_model = LangchainEmbedding(HuggingFaceEmbeddings(model_name='BAAI/bge-small-en')) # at this step the embedding files are downloaded
service_context = ServiceContext.from_defaults(embed_model=embed_model)



