#!/usr/bin/env python
# coding: utf-8

# In[332]:


wild_rows = df.copy()
wild_rows['mutant'] = wild_rows['wild_type']
wild_rows['mutant_1_letter'] = wild_rows['wild_type_1_letter']
wild_rows['delta_Tm_K'] = 0
wild_rows=wild_rows.drop_duplicates(subset='residue_number')
wild_rows


# In[333]:


df= pd.concat([df,wild_rows],ignore_index=True)
df


# In[334]:


# (Continue from previous code)

# Drop duplicates to avoid multiple mutations at the same site
wt_df = df[['chain', 'residue_number', 'wild_type_1_letter']].drop_duplicates()

# Sort by chain and residue number
wt_df = wt_df.sort_values(by=['chain', 'residue_number'])

# Group by chain and reconstruct the sequence
wt_sequences = wt_df.groupby('chain')['wild_type_1_letter'].apply(lambda x: ''.join(x))

# Print the sequence(s)
for chain, seq in wt_sequences.items():
    print(f"Wild-type sequence for chain {chain}:")
    print(seq)




# In[335]:


#seq as in testing set
our_model_masked_pos = pd.read_csv('5AZU_masked_position_predictions.csv')
ss=our_model_masked_pos['sequence'][0]
ss


# In[336]:


ss[20:]


# In[337]:


wt_df = df[['chain', 'residue_number', 'wild_type_1_letter']].drop_duplicates()

# Sort by chain and residue number
wt_df = wt_df.sort_values(by=['chain', 'residue_number'])

# Group by chain and reconstruct the sequence
wt_sequences = wt_df.groupby('chain')['wild_type_1_letter'].apply(lambda x: ''.join(x))

# Print the sequence(s)
for chain, seq in wt_sequences.items():
    print(f"Wild-type sequence for chain {chain}:")
    print(seq)



# In[338]:


df['residue_number']=df['residue_number']+20 #aligned to our seq


# In[339]:


df


# In[340]:


#ss[4:4+570]
misaligned=[i for i in range(len(seq)) if seq[i]!=ss[20:][i]]
misaligned


# In[341]:


df.to_csv('5AZU_HOTMusic_pred_cleaned.csv')


# In[342]:


our_model_masked_pos


# In[ ]:





# In[343]:


our_model_masked_pos.shape


# In[344]:


merged= pd.merge(our_model_masked_pos,df,left_on=['masked_position','predicted_aa','wt'],right_on=['residue_number','mutant_1_letter','wild_type_1_letter'],how='inner')
merged


# In[345]:


# merged= pd.merge(our_model_masked_pos,df,left_on=['masked_position','predicted_aa','wt'],right_on=['residue_number','mutant_1_letter','wild_type_1_letter'],how='left')
# merged


# In[346]:


merged.to_csv('5AZU_our_sugggested_aa_and_HOTMusic_pred_dtm.csv')


# In[ ]:





# In[ ]:





# In[ ]:




