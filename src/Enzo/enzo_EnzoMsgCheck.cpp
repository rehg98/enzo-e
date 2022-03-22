// See LICENSE_CELLO file for license and copyright information

/// @file     charm_EnzoMsgCheck.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-03-19
/// @brief    [\ref Charm] Declaration of the EnzoMsgCheck Charm++ message

#include "enzo.hpp"
#include "charm_enzo.hpp"

//----------------------------------------------------------------------

long EnzoMsgCheck::counter[CONFIG_NODE_SIZE] = {0};

//----------------------------------------------------------------------

EnzoMsgCheck::EnzoMsgCheck()
  : CMessage_EnzoMsgCheck(),
    is_local_(true),
    index_send_(),
    data_msg_(nullptr),
    buffer_(nullptr),
    block_name_(),
    block_lower_(),
    block_upper_(),
    block_size_(),
    tag_(),
    io_block_(),
    index_this_(),
    index_next_(),
    name_this_(),
    name_next_(),
    index_block_(),
    is_first_(),
    is_last_(),
    name_dir_()
{
  ++counter[cello::index_static()];
  cello::hex_string(tag_,TAG_LEN);
}

//----------------------------------------------------------------------

EnzoMsgCheck::~EnzoMsgCheck()
{
  --counter[cello::index_static()];
  CkFreeMsg (buffer_);
  buffer_=nullptr;
}

//----------------------------------------------------------------------

void EnzoMsgCheck::set_data_msg  (DataMsg * data_msg)
{
  if (data_msg_) {
    WARNING ("EnzoMsgCheck::set_data_msg()",
	     "overwriting existing data_msg_");
    delete data_msg_;
  }
  data_msg_ = data_msg;
}

//----------------------------------------------------------------------

void * EnzoMsgCheck::pack (EnzoMsgCheck * msg)
{
  // Return with buffer if already packed
  if (msg->buffer_ != nullptr) return msg->buffer_;

  int size = 0;

  // determine buffer size

  SIZE_OBJECT_TYPE(size,msg->index_send_);

  // data_msg_;
  int have_data = (msg->data_msg_ != nullptr);
  size += sizeof(int);
  if (have_data) {
    size += msg->data_msg_->data_size();
  }

  // Block name
  SIZE_STRING_TYPE(size,msg->block_name_);
  SIZE_ARRAY_TYPE(size,double,msg->block_lower_,3);
  SIZE_ARRAY_TYPE(size,double,msg->block_upper_,3);
  SIZE_ARRAY_TYPE(size,double,msg->block_upper_,3);

  SIZE_ARRAY_TYPE(size,char,msg->tag_,TAG_LEN+1);

  int have_io = (msg->io_block_ != nullptr);
  SIZE_SCALAR_TYPE(size,int,have_io);
  if (have_io) {
    SIZE_OBJECT_TYPE(size,*(msg->io_block_));
  }

  SIZE_OBJECT_TYPE(size,msg->index_this_);
  SIZE_OBJECT_TYPE(size,msg->index_next_);
  SIZE_STRING_TYPE(size,msg->name_this_);
  SIZE_STRING_TYPE(size,msg->name_next_);
  SIZE_SCALAR_TYPE(size,int,msg->index_block_);
  SIZE_SCALAR_TYPE(size,bool,msg->is_first_);
  SIZE_SCALAR_TYPE(size,bool,msg->is_last_);
  SIZE_STRING_TYPE(size,msg->name_dir_);

  //--------------------------------------------------

  // allocate buffer using CkAllocBuffer()

  char * buffer = (char *) CkAllocBuffer (msg,size);

  // serialize message data into buffer

  union {
    char * pc;
    int  * pi;
  };

  pc = buffer;

  SAVE_OBJECT_TYPE(pc,msg->index_send_);

  // data_msg_;
  have_data = (msg->data_msg_ != nullptr);
  (*pi++) = have_data;
  if (have_data) {
    pc = msg->data_msg_->save_data(pc);
  }

  // Block name
  SAVE_STRING_TYPE(pc,msg->block_name_);
  SAVE_ARRAY_TYPE(pc,double,msg->block_lower_,3);
  SAVE_ARRAY_TYPE(pc,double,msg->block_upper_,3);
  SAVE_ARRAY_TYPE(pc,double,msg->block_upper_,3);

  SAVE_ARRAY_TYPE(pc,char,msg->tag_,TAG_LEN+1);

  have_io = (msg->io_block_ != nullptr);
  SAVE_SCALAR_TYPE(pc,int,have_io);
  if (have_io) {
    SAVE_OBJECT_TYPE(pc,*(msg->io_block_));
  }

  SAVE_OBJECT_TYPE(pc,msg->index_this_);
  SAVE_OBJECT_TYPE(pc,msg->index_next_);
  SAVE_STRING_TYPE(pc,msg->name_this_);
  SAVE_STRING_TYPE(pc,msg->name_next_);
  SAVE_SCALAR_TYPE(pc,int,msg->index_block_);
  SAVE_SCALAR_TYPE(pc,bool,msg->is_first_);
  SAVE_SCALAR_TYPE(pc,bool,msg->is_last_);
  SAVE_STRING_TYPE(pc,msg->name_dir_);

  ASSERT2("EnzoMsgCheck::pack()",
	  "buffer size mismatch %ld allocated %d packed",
	  (pc - (char*)buffer),size,
	  (pc - (char*)buffer) == size);

  CkFreeMsg (msg);
  // Return the buffer
  return (void *) buffer;
}

//----------------------------------------------------------------------

EnzoMsgCheck * EnzoMsgCheck::unpack(void * buffer)
{

  // Allocate storage using CkAllocBuffer (not new!)

  EnzoMsgCheck * msg = (EnzoMsgCheck *) CkAllocBuffer (buffer,sizeof(EnzoMsgCheck));

  msg = new ((void*)msg) EnzoMsgCheck;

  msg->is_local_ = false;

  // de-serialize message data from input buffer into allocated message

  union {
    char   * pc;
    int    * pi;
  };

  pc = (char *) buffer;

  LOAD_OBJECT_TYPE(pc,msg->index_send_);

  // data_msg_
  int have_data = (*pi++);
  if (have_data) {
    msg->data_msg_ = new DataMsg;
    pc = msg->data_msg_->load_data(pc);
  } else {
    msg->data_msg_ = nullptr;
  }

  // Block name
  LOAD_STRING_TYPE(pc,msg->block_name_);
  LOAD_ARRAY_TYPE(pc,double,msg->block_lower_,3);
  LOAD_ARRAY_TYPE(pc,double,msg->block_upper_,3);
  LOAD_ARRAY_TYPE(pc,double,msg->block_upper_,3);

  LOAD_ARRAY_TYPE(pc,char,msg->tag_,TAG_LEN+1);

  int have_io;
  LOAD_SCALAR_TYPE(pc,int,have_io);
  if (have_io) {
    // create the correct IoBlock (IoBlock or IoEnzoBlock)
    msg->io_block_ = enzo::factory()->create_io_block();
    LOAD_OBJECT_TYPE(pc,*(msg->io_block_));
  }

  LOAD_OBJECT_TYPE(pc,msg->index_this_);
  LOAD_OBJECT_TYPE(pc,msg->index_next_);
  LOAD_STRING_TYPE(pc,msg->name_this_);
  LOAD_STRING_TYPE(pc,msg->name_next_);
  LOAD_SCALAR_TYPE(pc,int,msg->index_block_);
  LOAD_SCALAR_TYPE(pc,bool,msg->is_first_);
  LOAD_SCALAR_TYPE(pc,bool,msg->is_last_);
  LOAD_STRING_TYPE(pc,msg->name_dir_);

  // Save the input buffer for freeing later

  msg->buffer_ = buffer;

  return msg;
}

//----------------------------------------------------------------------

void EnzoMsgCheck::update (Data * data)
{
  // return if no data to update
  if (data_msg_ == nullptr) return;
  data_msg_->update(data,is_local_);

  if (!is_local_) {
    CkFreeMsg (buffer_);
    buffer_ = nullptr;
  }
}

//----------------------------------------------------------------------

Index EnzoMsgCheck::index_send()
{ return index_send_; }

//----------------------------------------------------------------------

void EnzoMsgCheck::set_index_send(Index index)
{ index_send_ = index; }

//----------------------------------------------------------------------

void EnzoMsgCheck::set_block (Block * block)
{
  block_name_ = block->name();
  block->data()->lower(block_lower_,block_lower_+1,block_lower_+2);
  block->data()->upper(block_upper_,block_upper_+1,block_upper_+2);
  block->data()->field().size(block_size_,block_size_+1,block_size_+2);

  io_block_ = enzo::factory()->create_io_block();
  io_block_->set_block(block);
}

//----------------------------------------------------------------------

void EnzoMsgCheck::del_block()
{
  delete data_msg_;
  data_msg_ = nullptr;
  delete io_block_;
  io_block_ = nullptr;
}

//----------------------------------------------------------------------

void EnzoMsgCheck::print (const char * msg)
{
  CkPrintf ("ENZO_MSG_CHECK====================\n");
  CkPrintf ("ENZO_MSG_CHECK tag %s %s\n",msg,tag_);
  CkPrintf ("ENZO_MSG_CHECK is_local %d\n",is_local_);
  int v3[3];
  index_send_.values(v3);
  CkPrintf ("ENZO_MSG_CHECK index_send values %d %d %d\n",v3[0],v3[1],v3[2]);
  CkPrintf ("ENZO_MSG_CHECK data_msg_ %p\n",(void *)data_msg_);
  CkPrintf ("ENZO_MSG_CHECK\n");
}
