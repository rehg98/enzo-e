// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodHello.cpp
/// @author   Ryan Golant (ryan.golant@columbia.edu)
/// @date     May 17, 2022
/// @brief    Implements the EnzoMethodHello class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodHello::EnzoMethodHello (double timeStep, int maxLevel)
  : Method(),
    timeStep_(timeStep),
    maxLevel_(maxLevel)
{
  //Refresh * refresh = cello::refresh(ir_post_);
  //cello::simulation()->refresh_set_name(ir_post_,name);
  
  ScalarDescr * scalar_descr_sync = cello::scalar_descr_sync();
  i_sync_restrict_ = scalar_descr_sync->new_value("hello:restrict");

  ScalarDescr * scalar_descr_void = cello::scalar_descr_void();
  for (int ic=0; ic<cello::num_children(); ic++) {
    i_msg_restrict_[ic] = scalar_descr_void->new_value("hello:msg_restrict");
  }

  n_ = 1;

  
}

//----------------------------------------------------------------------

void EnzoMethodHello::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | timeStep_;
  p | maxLevel_;
  p | i_sync_restrict_;
  p | n_;
  PUParray(p, i_msg_restrict_, 8);
}

//----------------------------------------------------------------------

void EnzoMethodHello::compute ( Block * block) throw()
{

  Sync * sync_restrict = psync_restrict(block);
  sync_restrict->set_stop(1 + cello::num_children());
  
  compute_ (block);
}



//----------------------------------------------------------------------

double EnzoMethodHello::timestep ( Block * block ) throw()
{
  // initialize_(block);

  return timeStep_;
}

//======================================================================

void EnzoMethodHello::compute_ (Block * block) throw()
{
  EnzoBlock * enzo_block = enzo::block(block);
  Sync * sync_restrict = psync_restrict(block);

  int level = enzo_block->level();
  //CkPrintf("level = %d, index = %d\n", level, block->thisIndex);

  //CkPrintf("first element = %d\n", i_msg_restrict_[0]);
  
  if (block->is_leaf()) {
    //CkPrintf("%d\n", maxLevel_);
    
    begin_cycle_ (enzo_block);
    
  } else {

    if (0 <= level && level < maxLevel_) {
      restrict_recv(enzo_block, nullptr);
    }
    

    else{
      CkPrintf("Hello!\n"); 
      }
    
  }
    
}


void EnzoMethodHello::begin_cycle_(EnzoBlock * enzo_block) throw()
{

  const int level = enzo_block->level();

  //CkPrintf("level = %d\n", level);

  if (level == 0) {
    CkPrintf("%d\n", n_); 
  }
  
  else {
    restrict (enzo_block);
  }

  enzo_block->compute_done();
}


void EnzoBlock::p_method_hello_restrict()
{
  EnzoMethodHello * method = 
    static_cast<EnzoMethodHello*> (this->method());

  method->restrict(this);
}


void EnzoMethodHello::restrict(EnzoBlock * enzo_block) throw()
{
 
  restrict_send (enzo_block);
}


void EnzoMethodHello::restrict_send(EnzoBlock * enzo_block) throw()
{

  
  FieldMsg * msg = pack_residual_(enzo_block);

  Index index_parent = enzo_block->index().index_parent(0);

  //CkPrintf("parent index = %d\n", index_parent);

  enzo::block_array()[index_parent].p_method_hello_restrict_recv(msg);
}


void EnzoBlock::p_method_hello_restrict_recv(FieldMsg * msg)
{
  EnzoMethodHello * method = 
    static_cast<EnzoMethodHello*> (this->method());

  method->restrict_recv(this,msg);
}


void EnzoMethodHello::restrict_recv
(EnzoBlock * enzo_block, FieldMsg * msg) throw()
{
  int level = enzo_block->level();
  //CkPrintf("In restrict_recv: level = %d, index = %d\n", level, enzo_block->thisIndex);
  //CkPrintf("start restrict_recv i_msg = %d\n", i_msg_restrict_[0]);  

  // Save field message from child
  if (msg != nullptr)
    {
      // print child_index(), print msg  
      //CkPrintf("Save msg: child_index = %d\n", msg->child_index());
      //CkPrintf("Message n: %d\n", msg->n);

      *pmsg_restrict(enzo_block,msg->child_index()) = msg;

      //CkPrintf("made it past pmsg_restrict!");
    }

  
  //CkPrintf("i_msg middle restrict_recv = %d\n", i_msg_restrict_[0]);
  
  // Continue if all expected messages received
  if (psync_restrict(enzo_block)->next() ) {

    //CkPrintf("all messages received\n");
    
    // Restore saved messages then clear
    for (int i=0; i<cello::num_children(); i++) {
      msg = *pmsg_restrict(enzo_block,i);
      *pmsg_restrict(enzo_block,i) = nullptr;
      // Unpack field from message then delete message
      unpack_residual_(enzo_block,msg);
    }
  
    begin_cycle_ (enzo_block);
  }

  //CkPrintf("end restrict_recv i_msg = %d\n", i_msg_restrict_[0]);

  
}


FieldMsg * EnzoMethodHello::pack_residual_(EnzoBlock * enzo_block) throw()
{
 
  // Create a FieldMsg for sending data to parent
 
  FieldMsg * msg  = new FieldMsg;
 
  msg->n = n_;
  char arr[1] = "";
  msg->a = arr;

  const int level = enzo_block->level();
  if (level > 0) 
  {
    enzo_block->index().child(level,msg->ic3,msg->ic3+1,msg->ic3+2);
  }

  //CkPrintf("%d, %d, %d\n", msg->ic3, msg->ic3+1, msg->ic3+2);

  return msg;

}


void EnzoMethodHello::unpack_residual_
(EnzoBlock * enzo_block,FieldMsg * msg) throw()
{
  // copy data from msg to this EnzoBlock

  n_ += 1;

  delete msg;
}
