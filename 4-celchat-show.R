rm(list = ls())
library(Seurat)
library(CellChat)
library(patchwork)
library(NMF)
library(ggalluvial)

load("1-cellchat.rda")
df.net <- subsetCommunication(cellchat)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

netAnalysis_signalingRole_scatter(cellchat)

table(cellchat@netP$pathways)
table(cellchat@idents)

netAnalysis_signalingRole_scatter(cellchat, signaling = "TGFb")

pathways.show <- cellchat@netP$pathways
pathways.show <- c("TGFb")
pathways.show <- c("VEGF")
netVisual_aggregate(cellchat,signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netAnalysis_contribution(cellchat, signaling = pathways.show)

plotGeneExpression(cellchat, signaling = c("VEGF"), enriched.only = T)
