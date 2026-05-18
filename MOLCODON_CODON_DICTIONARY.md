# MolCodon V2 - Codon Dictionary

Bu doküman, MolCodon V2 alfabesindeki tüm 3 harfli kodonların sınıflarını, anlamlarını ve neyi temsil ettiklerini özetler.

## 1. Yapısal Kontrol Kodonları (Structural Control)
Molekülün başlangıcını, bitişini ve halka birleşim noktalarını yönetir.
| Codon | Kategori | Anlamı / İşlevi |
|-------|----------|-----------------|
| `SCC` | Start | Molekül diziliminin başladığını belirtir. Daima ilk kodondur. |
| `SSS` | End | Molekül diziliminin bittiğini belirtir. Daima son kodondur. |
| `OCC` | Fusion | Kesişen (fused) halkalarda, atoma ikinci veya daha fazla kez gelindiğinde konur. |

## 2. Atom Kodonları (Elements)
Organik kimyadaki temel elementleri temsil eder.
| Codon | Element |
|-------|---------|
| `CCC` | Karbon (C) |
| `CCN` | Azot (N) |
| `CCO` | Oksijen (O) |
| `CCS` | Kükürt (S) |
| `CNC` | Flor (F) |
| `CNN` | Klor (Cl) |
| `CNO` | Brom (Br) |
| `CNS` | İyot (I) |
| `COC` | Fosfor (P) |
| `CON` | Bor (B) |

## 3. Bağ Kodonları (Bonds)
Atomlar arası bağlantının tipini belirtir.
| Codon | Bağ Tipi |
|-------|----------|
| `NCC` | Tekli (Single) |
| `NCN` | Çift (Double) |
| `NCO` | Üçlü (Triple) |
| `NCS` | Aromatik (Aromatic) |

## 4. Topoloji: Dallar (Branches)
Ana zincirden ayrılan ve geri dönen yan yolları belirtir. Aynı atoma bağlı birden fazla dal varsa farklı kodonlar (NNC, NOC vb.) sırayla kullanılarak dallar ayırt edilir.
| Açılış | Kapanış | İşlevi |
|--------|---------|--------|
| `NNC` | `NNN` | Dal 1 Açılış / Kapanış |
| `NOC` | `NON` | Dal 2 Açılış / Kapanış |
| `NOS` | `NOO` | Dal 3 Açılış / Kapanış |
| `NSC` | `NSN` | Dal 4 Açılış / Kapanış |

## 5. Topoloji: Halkalar (Rings)
Halkaların başladığı ve bittiği yerleri işaretler. Açılış ve kapanış kodonları eşleşerek (örn: NNO ile NNS) hangi halkanın kapandığını gösterir.
| Açılış | Kapanış |
|--------|---------|
| `NNO` | `NNS` |
| `NSO` | `NSS` |
| `OSO` | `OSS` |
| `OCN` | `ONO` |
| `OCO` | `ONS` |
| `OCS` | `OSC` |
| `ONC` | `OSN` |
| `ONN` | `SCN` |

## 6. Halka İçi Referans ve Konum Kodonları
Bir halka daha önce yazılmış başka bir halka ile kesişiyorsa (fused), eski halkanın indeksini ve atomun o halkadaki konumunu belirtmek için kullanılır.
- **Ring Reference (`SOC` - `SSO`):** Hangi halkaya bağlanıldığını gösterir (Halka 0, Halka 1 vb.).
- **Position (`COO` - `SNS`):** Atomun o halkada kaçıncı sırada olduğunu belirtir (0'dan 15'e kadar indeksler).

## 7. Atom Ek Açıklamaları (Atom Annotations)
Atomun fizikokimyasal veya üç boyutlu özelliklerini belirtmek için atom kodonunun hemen ardından yazılır.

### Yük (Charge)
| Codon | Anlamı |
|-------|--------|
| `CCX` | Nötr (0) |
| `CXN` | Pozitif (+1) |
| `CXS` | Yüksek Pozitif (+2 veya üstü) |
| `CXO` | Negatif (-1) |
| `CXX` | Yüksek Negatif (-2 veya altı) |

### Atom Stereokimyası
| Codon | Anlamı |
|-------|--------|
| `SXN` | R Konfigürasyonu |
| `SXO` | S Konfigürasyonu |

### Farmakofor Özellikleri
| Codon | Anlamı |
|-------|--------|
| `OXN` | Hidrojen Bağı Alıcı (HBA) |
| `OXC` | Hidrojen Bağı Verici (HBD) |

## 8. Bağ Ek Açıklamaları (Bond Annotations)
Bağın fizikokimyasal (esneklik) veya üç boyutlu özelliklerini belirtmek için bağ kodonunun hemen ardından yazılır.

### Esneklik / Hareketlilik (Mobility & Rotatability)
| Codon | Anlamı | Kapsamı |
|-------|--------|---------|
| `NXC` | Dönebilen (Rotatable) | Halka dışı, uçta olmayan tekli bağlar |
| `NCX` | Dönemeyen (Non-rotatable) | Halka dışı çift/üçlü bağlar veya uç (terminal) tekli bağlar |
| `NXS` | Halka Kısıtlamalı (Ring-constrained) | Aromatik olmayan halka içi bağlar |
| `NXO` | Aromatik Kilitli (Aromatic-locked) | Aromatik halka içi bağlar |

### Bağ Stereokimyası
| Codon | Anlamı |
|-------|--------|
| `SOX` | E (Trans) Konfigürasyonu |
| `SNX` | Z (Cis) Konfigürasyonu |

---

## Özet: Backbone (Omurga) Kavramı
"Backbone", kendi başına spesifik bir kodon türü değildir. Bir molekül encode edilirken;
- Dal (Branch) bloğu içine girmeyen,
- Halka (Ring) bloğu içine girmeyen,
ana ilerleyiş yolu üzerindeki **atomlar, bağlar ve bunların ek açıklamalarının** bütününe *Backbone* denir. BLAST benzerliği hesaplanırken bu omurga dizilimi Jaccard (N-gram) yöntemiyle karşılaştırılır.
