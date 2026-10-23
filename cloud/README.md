# Cloud run: SV calling + Exomiser with REMM

Runs `run_cloud.sh` on one AWS EC2 VM. It does five things:
1. Downloads the 85 GB of reads, the GRCh38 reference, the Exomiser data and REMM.
2. Aligns the reads.
3. Calls structural variants (Manta) and copy-number changes (CNVpytor).
4. Reruns Exomiser in whole-genome mode with REMM, once for SNVs/indels and once for the SV calls.
5. Packs the results into a tarball of a few MB.

**Estimated cost:** about 3.5 h × ~$4.30/h ≈ **$15**. **Your involvement:** about 15 minutes; the rest runs unattended.

Before starting, confirm that the hackathon data terms allow processing on a cloud VM you control. The VM and its disk are destroyed at step 6.

## 0. Check your CPU quota first (new accounts are often capped at 5–32 vCPUs)
AWS Console → **Service Quotas** → AWS services → **Amazon EC2** → search *Running On-Demand Standard (A, C, D, H, I, M, R, T, Z) instances*. Make sure the **us-east-1** region is selected.
- **96 or more:** go ahead.
- **32–95:** use `m7i.8xlarge` below instead. It's ~2× slower and about the same cost.
- **Under 32:** click *Request increase at account level* and ask for **96**. Approval usually takes minutes to a day, so do this today.

## 1. Launch the VM (AWS Console → EC2 → Launch instance)
| Setting | Value |
|---|---|
| Region | **us-east-1** (fastest route to Hugging Face) |
| Name | `mva-sv` |
| AMI | **Ubuntu Server 24.04 LTS**, 64-bit (x86) |
| Instance type | **c7i.24xlarge** (96 vCPU, 192 GiB). If AWS refuses it, m7i.8xlarge (32 vCPU) also works, just ~2× slower |
| Key pair | *Create new key pair* → name `mva-sv`, RSA, `.pem` → it downloads |
| Network | Allow SSH traffic from **My IP** only |
| Storage | **600 GiB gp3**, IOPS 6000, throughput 500 MB/s, *Delete on termination* = Yes |

Click *Launch*, then copy the instance's **Public IPv4 address**.

## 2. Copy the script up (from this repo, on your Mac)
```bash
mv ~/Downloads/mva-sv.pem ~/.ssh/ && chmod 400 ~/.ssh/mva-sv.pem
IP=<public-ip>
scp -i ~/.ssh/mva-sv.pem cloud/run_cloud.sh ubuntu@$IP:
```

## 3. Start it
```bash
ssh -i ~/.ssh/mva-sv.pem ubuntu@$IP
read -rs HF_TOKEN && export HF_TOKEN     # paste your Hugging Face token, then Enter (nothing is shown)
nohup bash run_cloud.sh > run.log 2>&1 &
exit
```
The `nohup` keeps the job running if your connection drops, so you can close the laptop.

## 4. Check progress (any time)
```bash
ssh -i ~/.ssh/mva-sv.pem ubuntu@$IP 'tail -5 run.log; ls data/done'
```
The steps run in this order: `tools`, `downloads`, `exo_unpack`, `align`, `manta`, `cnv`, `exomiser`, `summary`. `align` is the long one, at roughly 1–1.5 h. REMM downloads in the background, and only `exomiser` waits for it.

**If it fails:** fix the problem and rerun the same `nohup` command. Finished steps are skipped.

## 5. Bring the results home
```bash
scp -i ~/.ssh/mva-sv.pem ubuntu@$IP:mva_results.tar.gz results/ && tar -xzf results/mva_results.tar.gz -C results/ && mv results/results results/cloud
```
Then ask Claude to review `results/cloud/SUMMARY.txt`.

## 6. Destroy the VM (this deletes all data on it)
EC2 → Instances → `mva-sv` → *Instance state* → **Terminate**. Then delete the `mva-sv` key pair under EC2 → Key Pairs.
